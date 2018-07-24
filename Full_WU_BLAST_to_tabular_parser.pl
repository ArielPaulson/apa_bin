#!/usr/bin/perl
use strict;

## WU-blastn output, e.g.:
## pressdb MM48_genome.fa
## wublastn -d MM48_genome.fa -i V1RD14_consensus1.txt -m 8 -W 3 > all_genome_hits_wu_W.txt

my (%data, %adata, %subjects, %returns, %posoffset, %laboffset, @chars);
my ($line, $alignline, $alignflag, $offset, $acount, $qflag, $sflag, $query, $subject, $file, $qcount);
my %strands = ('Plus' => '+', 'Minus' => '-');
my @QSFields1 = qw/ SCORE E_VAL IDTPCT LENGTH QLPCT QSTART QEND QSTRAND /;
my @QSFields2 = qw/ SSTART SEND SSTRAND /;

if (@ARGV) {
	$file = $ARGV[0];
	die "$file not found!" unless (-e $file);
} else {
	print "Enter a full BLAST output file to parse:\n";
	{
		chomp($file = <STDIN>);
		unless (-e $file) {
			print "$file not found!  Retry:\n";
			redo;
		}
	}
}

#$posoffset{$query} = 0;
open IN, $file;
while (my $line = <IN>) {
	$_ =~ s/[\n\r]//g;
	if ($line =~ /^Query=\s+(\S+)/) {
		$qcount++;
		$query = $1;
		$returns{QUERY}{$query} = $qcount;		# count returning queries
		$laboffset{$query} = length($query);	# seed the pileup offset maximum value
		$qflag = 1;			# next line has the query length
		$acount = 0;			# reset alignment counter
	} elsif ($qflag == 1) {
		$qflag = 0;
		$data{$query}{LENGTH} = $1 if ($line =~ /^\s+\(([,\d]+) letters\)/);
		$data{$query}{LENGTH} =~ s/,//g;
	} elsif ($sflag == 1) {
		$sflag = 0;
		$data{$subject}{LENGTH} = $1 if ($line =~ /^\s+Length = ([,\d]+)/);
		$data{$subject}{LENGTH} =~ s/,//g;
	} elsif ($line =~ /^>(.*)/) {
		($subject = $1) =~ s/\s+$//;		# clip trailing whitespace, if any
		$returns{SUBJECT}{$subject}++;	# count hits to subjects
		$sflag = 1;			# next line has the subject length
		$alignflag = 0;			# have not yet seen alignments for this subject
		$laboffset{$query} = length($subject) if (length($subject) > $laboffset{$query});  # retain maximum value for offset (if alignment will be shown)
	} elsif ($line =~ /^\s+Length =\s+(\d+)/) {
		$data{$subject}{LENGTH} = $1;
		$sflag = 0;
	} elsif ($line =~ /^\s+Score = (\d+) \(([\.\d]+) bits\), Expect = ([-e\.\d]+), P = ([-e\.\d]+)/) {
		$acount++;			# every alignment always starts with a "Score" line
		$returns{QHITS}{$1}++;		# count hits to queries
		$adata{$query}{$acount}{SCORE} = $1;
		$adata{$query}{$acount}{E_VAL} = $3;
		$adata{$query}{$acount}{P_VAL} = $4;
		$adata{$query}{$acount}{SBJCT} = $subject;
		if ($1 > $data{$query}{TOPSCORE}) {
			$data{$query}{TOPSCORE} = $1; 
			$data{$query}{TOPHIT} = $acount;
		}
		# prep for later
		$adata{$query}{$acount}{QSTART} = 9E9;
		$adata{$query}{$acount}{SSTART} = 9E9;
	} elsif ($line =~ /^\s+Identities = (\d+)\/(\d+) \((\d+)%\), Positives = (\d+)\/(\d+) \((\d+)%\), Strand = (\w+) \/ (\w+)/) {
		$adata{$query}{$acount}{QIDENT} = $1;
		$adata{$query}{$acount}{SIDENT} = $1;
		$adata{$query}{$acount}{LENGTH} = $2;
		$adata{$query}{$acount}{IDTPCT} = $3;
		$adata{$query}{$acount}{QLPCT} = sprintf("%0.0f", 100*$3/$data{$query}{LENGTH});
		$adata{$query}{$acount}{SLPCT} = sprintf("%0.0f", 100*$3/$data{$subject}{LENGTH});
		# ignoring $4,$5,$6 -- these are normally redundant to $1,$2,$3
		$adata{$query}{$acount}{QSTRAND} = $strands{$7};
		$adata{$query}{$acount}{SSTRAND} = $strands{$8};
	} elsif ($line =~ /^Query:\s+(\d+)(\s+)([-\w]+)\s+(\d+)/) {
		$adata{$query}{$acount}{QSTART} = $1 if ($1 < $adata{$query}{$acount}{QSTART});
		$adata{$query}{$acount}{QEND} = $4 if ($4 > $adata{$query}{$acount}{QEND});
#		$adata{$query}{$acount}{QSEQ} .= $3;
#		$posoffset{$query} = length($1) if (length($1) > $posoffset{$query});
#		push @{ $adata{$query}{$acount}{ALIGN} }, "$line\n";
#		$offset = 7 + length($1) + length($2);	# character position where alignment symbols begin (may start w/ space)
#		$alignflag = 1;				# have seen at least one alignment for this subject
		$alignline = 1;				# next line has the alignment symbols
	} elsif ($alignline == 1) {
		$alignline = 0;
#		my @chars = split //, $line;
#		$adata{$query}{$acount}{IDALN} .= join "", @chars[$offset..$#chars];
#		push @{ $adata{$query}{$acount}{ALIGN} }, "$line\n";
	} elsif ($line =~ /^Sbjct:\s+(\d+) ([-\w]+) (\d+)/) {
		$adata{$query}{$acount}{SSTART} = $1 if ($1 < $adata{$query}{$acount}{SSTART});
		$adata{$query}{$acount}{SEND} = $3 if ($3 > $adata{$query}{$acount}{SEND});
#		$adata{$query}{$acount}{SSEQ} .= $2;
#		$posoffset{$query} = length($1) if (length($1) > $posoffset{$query});
#		push @{ $adata{$query}{$acount}{ALIGN} }, "$line\n\n";
	}
}
close IN;

open OUT, "> tabular_parsed_$file";
print OUT join "\t", (@QSFields1, 'SBJCT', @QSFields2), "\n";
foreach my $query (sort { $returns{QUERY}{$a} <=> $returns{QUERY}{$b} } keys %{ $returns{QUERY} }) {
	foreach my $align (sort {$a <=> $b} keys %{ $adata{$query} }) {
		my $line1 = join "\t", map { $adata{$query}{$align}{$_} } @QSFields1;
		my $line2 = join "\t", map { $adata{$query}{$align}{$_} } @QSFields2;
		print OUT "$line1\t$adata{$query}{$align}{SBJCT}\t$line2\n";
	}
}
close OUT;
exit;

