#!/usr/bin/env perl
use Getopt::Long;

my ($bfile, $qfile, $sfile, $ofile, $lfile, $minmatch, $minidt, $gapless, $best);
GetOptions(
    "b=s" => \$bfile, 
    "q=s" => \$qfile, 
    "s=s" => \$sfile, 
    "o=s" => \$ofile, 
    "l=s" => \$lfile, 
    "mm=s" => \$minmatch, 
    "mi=s" => \$minidt, 
    "gapless" => \$gapless,
    "best" => \$best
);

my ($monoQ, $monoS, $numQ, $numS);
if ($qfile) {
    if ($qfile =~ /^\d+$/) {  # single-value only
	$monoQ = $sfile;
    } else {
	open IN, $qfile;
	while (<IN>) {
	    chomp;
	    my ($query, $length) = split /\t/, $_;
	    $hdata{Q}{$query} = $length;
	}
	close IN;
    }
    $numQ = scalar keys %{ $hdata{Q} };
}

if ($sfile) {
    if ($sfile =~ /^\d+$/) {
	$monoS = $sfile;
    } else {
	open IN, $sfile;
	while (<IN>) {
	    chomp;
	    my ($subject, $length) = split /\t/, $_;
	    $hdata{S}{$subject} = $length;
	}
	close IN;
    }
    $numS = scalar keys %{ $hdata{S} };
}

my (%inQ, %inS, %hitQ, %hitS);

my %scored;
open IN3, $bfile;
while (<IN3>) {
    chomp($line = $_);
    if ($. < 3) {
        next;
    } elsif ($. < 6) {
        push @header, $line;
        next;
    }
    my @data = split /\t/, $line;
    my ($identbp, $misses, $reps, $Ns, $Qgaps, $Qgapbp, $Tgaps, $Tgapbp, $strand, $query, $qlen, $qpos1, $qpos2, $subject, $slen, $spos1, $spos2, $Nblocks, $blocksizes, $qstarts, $tstarts) = @data;
    $inQ{$query} = 1;
    $inS{$subject} = 1;
    next if ($gapless && ($Qgaps || $Tgaps));
    next if ($minmatch && $identbp < $minmatch);
    next if ($minidt && $identbp / $mlen < $minidt);
    my $lenQ = $hdata{Q}{$query} || $monoQ;
    my $lenS = $hdata{S}{$subject} || $monoS;
    my ($mrangeQ, $mrangeS) = (abs($qpos2 - $qpos1) + 1, abs($spos2 - $spos1) + 1);
    my $mlen = $mrangeQ < $mrangeS ? $mrangeQ : $mrangeS;
    my $identpct = sprintf("%0.3f", $identbp / $mlen);	                # clean up before printing
    my $Ngaps = $Qgaps + $Sgaps;
    my $Ngapbp = $Qgapbp + $Sgapbp;
    my ($covgQ, $covgS);
    $covgQ = sprintf("%0.3f", $mrangeQ/$lenQ) if $lenQ;	            # clean up before printing
    $covgS = sprintf("%0.3f", $mrangeS/$lenS) if $lenS;	            # clean up before printing
    $newline = join "\t", ($line, $identpct, $Ngaps, $Ngapbp, $mrangeQ, $lenQ, $covgQ, $mrangeS, $lenS, $covgS);
    $linecount++;
    if ($best) {
	$scored{$query}{$identbp}{$Qgapbp}{$newline} = $subject;
    } else {
	push @output, "$newline\n";
	$hitQ{$query} = 1;
        $hitS{$subject} = 1;
    }
}
close IN;

$numQ = scalar keys %inQ unless $numQ;
$numS = scalar keys %inS unless $numS;

if ($best) {
    foreach my $query (sort keys %scored) {
	$hitQ{$query} = 1;
	my $besthits;
	my $topscore = [sort {$b <=> $a} keys %{ $scored{$query} }]->[0];
	my $mingap = [sort {$a <=> $b} keys %{ $scored{$query}{$topscore} }]->[0];
	foreach my $line (keys %{ $scored{$query}{$topscore}{$mingap} }) {
	    $subject = $scored{$query}{$topscore}{$mingap}{$line};
	    $hitS{$subject} = 1;
	    push @output, "$line\n";
	    $allbest++;
	}
    }
}

my $inQ = scalar keys %inQ;
my $inS = scalar keys %inS;
my $hitQ = scalar keys %hitQ;
my $hitS = scalar keys %hitS;

my $lostQ = 0;  # guarantee printable 
if ($qfile) {
    my @lostQ;
    foreach my $query (keys %{ $hdata{Q} }) {
        push @lostQ, $query unless $hitQ{$query};
    }
    $lostQ = scalar @lostQ;
    my $lname;
    if ($lfile =~ /^[yY]$/) {
	$lname = "lost_${qfile}_$bfile.txt";
    } elsif ($lfile) {
	$lname = $lfile;
    }
    if ($lname) {
	open OUT2, "> $lname";
	print OUT2 "$_\n" foreach @lostQ;
	close OUT2;
    }
}

unless ($ofile) {
    my ($bname) = ($bfile =~ /([^\/]+)$/);
    $ofile = "best_alignments_$bname";
}
open OUT1, "> $ofile";
print OUT1 "$header[0]\tIdentPct\tQ_Range\tS_Range\tQ_Cov\tS_Cov\n$header[1]\n$header[2]\n";
print OUT1 @output;
close OUT1;
print "$numQ queries | $inQ blatted | $hitQ passed | $lostQ lost\n$numS subjects | $inS blatted | $hitS passed\n$linecount total alignments | $allbest best\n";
print "Process complete!\n";
exit;


