#!/usr/bin/perl
use strict;

my (%data, %adata, %subjects, %returns, %posoffset, %laboffset, @chars);
my ($line, $alignline, $alignflag, $offset, $acount, $qflag, $query, $subject, $file, $qcount);
my %strands = ('Plus' => '+', 'Minus' => '-');
my @QSFields1 = qw/ SCORE E_VAL IDTPCT GAPPCT LENGTH QLPCT SLPCT QIDENT QGAPS SIDENT SGAPS /;
my @QSFields2 = qw/ QSTART QEND QSTRAND SSTART SEND SSTRAND /;

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
		print "$1\n";
		$returns{QUERY}{$query} = $qcount;		# count returning queries
		$laboffset{$query} = length($query);	# seed the pileup offset maximum value
		$qflag = 1;			# next line has the query length
		$acount = 0;			# reset alignment counter
	} elsif ($qflag == 1) {
		$qflag = 0;
		$data{$query}{LENGTH} = $1 if ($line =~ /^\s+\((\d+) letters\)/);
	} elsif ($line =~ /^>(.*)/) {
		($subject = $1) =~ s/\s+$//;		# clip trailing whitespace, if any
		$returns{SUBJECT}{$subject}++;	# count hits to subjects
		$alignflag = 0;			# have not yet seen alignments for this subject
		$laboffset{$query} = length($subject) if (length($subject) > $laboffset{$query});  # retain maximum value for offset (if alignment will be shown)
	} elsif ($line =~ /^\s+Length =\s+(\d+)/) {
		$data{$subject}{LENGTH} = $1;
	} elsif ($line =~ /^\s+Score =\s+([\.\d]+) bits \((\d+)\), Expect =\s+([-e\.\d]+)/) {
		$acount++;			# every alignment always starts with a "Score" line
		$returns{QHITS}{$1}++;		# count hits to queries
		$adata{$query}{$acount}{SCORE} = $1;
		$adata{$query}{$acount}{LENGTH} = $2;
		$adata{$query}{$acount}{E_VAL} = $3;
		$adata{$query}{$acount}{SBJCT} = $subject;
		$adata{$query}{$acount}{QLPCT} = sprintf("%0.0f", 100*$2/$data{$query}{LENGTH});
		$adata{$query}{$acount}{SLPCT} = sprintf("%0.0f", 100*$2/$data{$subject}{LENGTH});
		if ($1 > $data{$query}{TOPSCORE}) {
			$data{$query}{TOPSCORE} = $1; 
			$data{$query}{TOPHIT} = $acount;
		}
		# prep for later
		$adata{$query}{$acount}{QSTART} = 9E9;
		$adata{$query}{$acount}{SSTART} = 9E9;
	} elsif ($line =~ /^\s+Identities =\s+(\d+)\/(\d+) \(([\.\d]+)%\)/) {
		$adata{$query}{$acount}{QIDENT} = $1;
		$adata{$query}{$acount}{SIDENT} = $2;
		$adata{$query}{$acount}{IDTPCT} = $3;
		if ($line =~ /Gaps =\s+(\d+)\/(\d+) \(([\.\d]+)%\)/) {
			$adata{$query}{$acount}{QGAPS} = $1;
			$adata{$query}{$acount}{SGAPS} = $2;
			$adata{$query}{$acount}{GAPPCT} = $3;
		} else {
			$adata{$query}{$acount}{QGAPS} = 0;
			$adata{$query}{$acount}{SGAPS} = 0;
			$adata{$query}{$acount}{GAPPCT} = 0;
		}
	} elsif ($line =~ /^\s+Strand =\s+(\w+) \/ (\w+)/) {
		$adata{$query}{$acount}{QSTRAND} = $strands{$1};
		$adata{$query}{$acount}{SSTRAND} = $strands{$2};
	} elsif ($line =~ /^Query: (\d+)(\s+)([-\w]+)\s+(\d+)/) {
		$adata{$query}{$acount}{QSTART} = $1 if ($1 < $adata{$query}{$acount}{QSTART});
		$adata{$query}{$acount}{QEND} = $4 if ($4 > $adata{$query}{$acount}{QEND});
		$adata{$query}{$acount}{QSEQ} .= $3;
#		$posoffset{$query} = length($1) if (length($1) > $posoffset{$query});
#		push @{ $adata{$query}{$acount}{ALIGN} }, "$line\n";
#		$offset = 7 + length($1) + length($2);	# character position where alignment symbols begin (may start w/ space)
#		$alignflag = 1;				# have seen at least one alignment for this subject
		$alignline = 1;				# next line has the alignment symbols
	} elsif ($alignline == 1) {
#		$alignline = 0;
#		my @chars = split //, $line;
#		$adata{$query}{$acount}{IDALN} .= join "", @chars[$offset..$#chars];
#		push @{ $adata{$query}{$acount}{ALIGN} }, "$line\n";
	} elsif ($line =~ /^Sbjct: (\d+)\s+([-\w]+)\s+(\d+)/) {
		$adata{$query}{$acount}{SSTART} = $1 if ($1 < $adata{$query}{$acount}{SSTART});
		$adata{$query}{$acount}{SEND} = $3 if ($3 > $adata{$query}{$acount}{SEND});
		$adata{$query}{$acount}{SSEQ} .= $2;
#		$posoffset{$query} = length($1) if (length($1) > $posoffset{$query});
#		push @{ $adata{$query}{$acount}{ALIGN} }, "$line\n\n";
	}
}
close IN;

open OUT, "> tabular_parsed_$file";
print OUT join "\t", (@QSFields1, @QSFields2), "\n";
foreach my $query (sort { $returns{QUERY}{$a} <=> $returns{QUERY}{$b} } keys %{ $returns{QUERY} }) {
	print "$query\n";
	foreach my $align (sort {$a <=> $b} keys %{ $adata{$query} }) {
		print "\t$align\n";
		my $line = join "\t", map { $adata{$query}{$align}{$_} } (@QSFields1, @QSFields2);
		print OUT "$line\n";
	}
}
close OUT;
exit;

