#!/usr/bin/env perl
use strict;

# Intended for analyzing a bowtie2-aligned set of trimmed reads.
# Specifically, reads aligned to one of apa's housekeeping-ncRNA references.
# Outputs a 3-column dataset for analysis in R: read length, ncRNA class (or '*' for unaligned), and 1/0 for multiread status.
# R analysis reports totals, alignability and multiplicity stats per read length, and read-length vs alignability/multireadness per ncRNA class.

my $bam = $ARGV[0];
my %hash;

if ($bam =~ /\.bam$/i) {
    open IN, "samtools view $bam |" or die "$0: Cannot 'samtools view' BAM file '$bam': $!\n";
} elsif ($bam =~ /\.sam$/i) {
    open IN, $bam or die "$0: Cannot open SAM file '$bam': $!\n";
} else {
    die "$0: Input file must have extension 'sam' or 'bam', not case sensitive.\n";
}
print "ReadLen\tClass\tMulti\n";
while (<IN>) {
    next if $_ =~ /^@/;
    my @data = split /\t/, $_, 11;
    my $mapq = $data[4];  # primary alignment score
    my ($sas) = ($data[10] =~ /XS:i:(\d+)/);  # best secondary-alignment score
    my $multi = $sas >= $mapq-1 ? 1 : 0;  # if secondary score == primary score (or primary-1), assert that there were > 1 equally-good alignments
    my $class = (split /\|/, $data[2])[0];  # type of ncRNA (or '*' if unaligned)
    my $readlen = length($data[9]);  # trimmed read length
    print "$readlen\t$class\t$multi\n";
}
close IN;
exit;
