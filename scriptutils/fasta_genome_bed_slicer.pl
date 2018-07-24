#!/n/site/inst/Linux-x86_64/sys/bin/perl

if (@ARGV) {
	die "Bad command line!\n" if scalar @ARGV != 2;
	($fadir, $bed) = @ARGV;
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
	
	print "\nEnter a bed file with locations:\n";
	{
		chomp($bed = <STDIN>);
		unless (-e $bed) {
			print "$bed does not exist!  Try again:\n";
			redo;
		}
	}
}

open IN, $bed;
while (<IN>) {
	$_ =~ s/[\n\r]//g;
    next if $_ =~ /^track/;
	next unless $_;		# no blank lines
	my ($chr, $start, $end, $id, $score, $strand) = split "\t", $_;
	my ($length1, $length2);
    $id = $. unless $id;
    $strand = '+' unless $strand;
	$allids{$id} = 0;	# gets 1 later if extracted successfully
	$groups{$chr}{$id} = [$start, $end, $strand];
}
close IN;

open OUT, "> $bed.fa";
foreach my $chr (sort keys %groups) {
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

open IDS, "> extraction_success_$bed";
print IDS "ID\tSuccess\n";
print IDS "$_\t$allids{$_}\n" foreach keys %allids;
close IDS;
exit;

sub extract {
	my ($chr, $seqref) = @_;
	foreach my $id (keys %{ $groups{$chr} }) {
		my ($start, $end, $strand) = @{ $groups{$chr}{$id} };
		my $seq;
		if ($start <= $end) {
			$seq = substr($$seqref, $start-1, $end-$start+1);
		} elsif ($start > $end) {	# wraps around breakpoint of circular chr, e.g. MT or plasmid
			my $seqA = substr($$seqref, $end-1, length($$seqref)-$end+1);
			my $seqB = substr($$seqref, 0, $start);
			$seq = $seqA.$seqB;
		}
		if ($seq) {
			$seq = &revcomp($seq) if ($strand eq '-' | $strand == -1);
			$seq = &blockify($seq);
			$allids{$id} = 1;
			print OUT ">$id\n$seq";
		}
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
