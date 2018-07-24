#!/n/site/inst/Linux-x86_64/sys/bin/perl

if (@ARGV) {
	die "Bad command line!\n" if scalar @ARGV != 2;
	($fadir, $psl) = @ARGV;
} else {
	print "Enter a path to chromosome fastas: ('.' if here)\n";
	{
		chomp($fadir = <STDIN>);
		unless (-d $fadir) {
			print "$fadir does not exist!  Try again:\n";
			redo;
		}
		$fadir = undef if $fadir eq '.';
	}
	
	print "\nEnter a psl file:\n";
	{
		chomp($psl = <STDIN>);
		unless (-e $psl) {
			print "$psl does not exist!  Try again:\n";
			redo;
		}
	}
}

my $active;
open IN, $psl;
while (<IN>) {
    $_ =~ s/[\n\r]+$//;
	$active = $. if $_ !~ /[^-\s]/;  # if fails to bind and non-hyphen, nonspace, then must be last header line
    next unless $. > $active;
	$in++;
    my @data = split /\t/, $_;
    my ($match, $misses, $repmatch, $Ns, $qgaps, $qgapbp, $sgaps, $sgapbp, $strand, $name, $length, $qpos1, $qpos2, $chr, $chrlen, $spos1, $spos2, $blocks, $blocksizes, $blockstarts, $chrstarts) = split /\t/, $_;
	$allids{$name} = 0;	# gets 1 later if extracted successfully
    $genestrand{$name} = $strand;
    my ($start, $end) = (sort {$a <=> $b} ($spos1, $spos2));
    if ($chr =~ /_random/) {
		(my $chr1 = $chr) =~ s/_random//;
		$chr1 =~ /\D/ ? ($alpha{$chr} = 1) : ($numeric{$chr} = 1);
    } else {
		$chr =~ /\D/ ? ($alpha{$chr} = 1) : ($numeric{$chr} = 1);
    }
    $blocksizes =~ s/,$//;
    $chrstarts =~ s/,$//;
    my @blocksizes = split /,/, $blocksizes;
    my @chrstarts = split /,/, $chrstarts;
    foreach my $i (0..$#blocksizes) {
		push @{ $groups{$chr}{$name} }, [$chrstarts[$i]+1, $chrstarts[$i]+$blocksizes[$i]];
	}
}
close IN;

my @chrs = ( (sort {$a <=> $b} keys %numeric), (sort keys %alpha) );

open OUT, "> $psl.fa";
foreach my $chr (@chrs) {
	($fadir) ? ($fa1 = "$fadir/$chr.fa") : ($fa1 = "$chr.fa"); 
	($fadir) ? ($fa2 = "$fadir/$chr.fasta") : ($fa2 = "$chr.fasta"); 
	my $nseq = scalar keys %{ $groups{$chr} };
	if (-e $fa1) {
		open FA, $fa1;
	} elsif (-e $fa2) {
		open FA, $fa2;
	} else {
		print "No obvious fasta files for chromosome $chr found!  Skipping $nseq entries...\n";
		next;
	}
	
	print "Reading $chr.fa...\n";
	my $chrseq;
	while (<FA>) {
		$_ =~ s/[\n\r]//g;
		$chrseq .= $_ if ($_ && $_ !~ /^>/);	# skip headers and the occasional blank line
	}
	close FA;
	print length($chrseq), " bases read.\n";

	print "Extracting $nseq sequences...\n";
	&extract($chr, \$chrseq);
}
close OUT;

open IDS, "> extraction_success_$psl";
print IDS "ID\tSuccess\n";
print IDS "$_\t$allids{$_}\n" foreach keys %allids;
close IDS;
exit;

sub extract {
	my ($chr, $seqref) = @_;
	foreach my $gene (keys %{ $groups{$chr} }) {
        my $geneseq;
		foreach my $exon (@{ $groups{$chr}{$gene} }) {  
            my ($start, $end) = @$exon;
            my $seq;
            if ($start <= $end) {   # normal
                $seq = substr($$seqref, $start-1, $end-$start+1);
#            } elsif ($start > $end) {	# wraps around breakpoint of circular chr, e.g. MT or plasmid
#                my $seqA = substr($$seqref, $end-1, length($$seqref)-$end+1);
#                my $seqB = substr($$seqref, 0, $start);
#                $seq = $seqA.$seqB;
            }
            if ($seq) {
                $seq = &revcomp($seq) if $genestrand{$gene} eq '-';
                $seq = &blockify($seq);
                $geneseq .= $seq;
                $allids{$gene} = 1;
            } else {
                $geneseq .= '*';  # lost sequence marker
            }
        }
        print OUT ">$gene\n$geneseq";
	}
}

sub revcomp {
	my $SEQ = shift;
	$SEQ =~ tr/ACGT/TGCA/;
	return reverse($SEQ);
}

sub blockify {
	my $SEQ = shift;
	my ($seqblock, $start);
	my $blocks = length($SEQ) / 50;
	$blocks++ if (length($SEQ) % 50 != 0);
	for (my $i = 1; $i <= $blocks; $i++) {
		$seqblock .= substr($SEQ, $start, 50)."\n";
		$start += 50;
	}
	return $seqblock;
}
