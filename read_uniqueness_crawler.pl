#!/usr/bin/perl

my $outfile = "all.read.uniqueness.txt";
my $temp = "temp.$$";

foreach my $dir (glob "*") {
    next unless -d $dir;
    chdir $dir;
    foreach my $bam (glob "*.bam") {
	my $tag = "$dir/$bam";
	push @tags, $tag;
	chomp(my $now = `date`);
	print "Processing $tag: $now\n";
	system "samtools view $bam | cut -f1,10 > $temp";
	my (%reads, %seqs);
	open IN, $temp;
	while (<IN>) {
	    chomp;
	    my ($read, $seq) = split /\t/, $_, 2;
	    $stats{$tag}{ALIGNS}++;
	    $reads{$read}++;
	    $seqs{$seq}++;
	}
	close IN;
	unlink $temp;
	$stats{$tag}{READS} = scalar keys %reads;
	$stats{$tag}{SEQS} = scalar keys %seqs;
    }
    chdir '..';
}
open OUT, "> $outfile" or die "Couldn't write output to '$outfile': $!\n";
print OUT "\t$_" foreach @tags;
print OUT "\n";
foreach my $type (qw/ ALIGNS READS SEQS /) {
    print OUT $type;
    print OUT "\t$stats{$_}{$type}" foreach @tags;
    print OUT "\n";
}
close OUT;
exit;

