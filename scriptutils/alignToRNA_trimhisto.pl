#!/usr/bin/env perl
use strict;

# Intended for analyzing a bowtie2-aligned set of trimmed reads.
# Specifically, reads aligned to one of apa's housekeeping-ncRNA references.
# Outputs a 3-column dataset for analysis in R: read length, ncRNA class (or '*' for unaligned), and 1/0 for multiread status.
# R analysis reports totals, alignability and multiplicity stats per read length, and read-length vs alignability/multireadness per ncRNA class.

my ($samp, $outfile) = @ARGV;
my $hskpbam = "$samp/$samp.bam";  # ncRNA bam
my $genobam = "$samp/$samp.genome/accepted_hits.bam";  # full tophat genome bam, NOT primary

my %data;

my %classes = ('*',0, 'miRNA',1, 'misc_RNA',2, 'Mt_rRNA',3, 'Mt_tRNA',4, 'rRNA',5, 'snoRNA',6, 'snRNA',7, 'tRNA',8);

open my $OUT, '>', $outfile or die "Cannot write to '$outfile': $!\n";

if (-e $genobam) {
    
    ## process Hskp_ncRNA bam
    ## single alignment per read
    open my $IN2, '-|', "samtools view $hskpbam" or die "$0: Cannot 'samtools view' BAM file '$hskpbam': $!\n";
    while (<$IN2>) {
	my @fields = split /\t/, $_, 11;
	$data{$fields[0]}{C} = $classes{ (split /\|/, $fields[2])[0] };  # type of ncRNA (or '*' if unaligned), as integer 0-8
    }
    close $IN2;
    
    ## process Genome bam
    ## multiple alignments per read
    open my $IN1, '-|', "samtools view $genobam" or die "$0: Cannot 'samtools view' BAM file '$genobam': $!\n";
    while (<$IN1>) {
	my @fields = split /\t/, $_, 11;
	$data{$fields[0]}{L} = length($fields[9]);  # trimmed read length
	$data{$fields[0]}{I} = $. unless $data{$fields[0]}{I};
	$data{$fields[0]}{N}++;
    }
    close $IN1;
    
    ## output
#    print "Read\tTrimLen\tClass\tAligns\n";
    print "TrimLen\tClass\tAligns\n";
    foreach my $idx (sort {$a <=> $b} keys %data) {
#	print $OUT join("\t", map {$data{$idx}{$_}} qw/ I L C N /), "\n";
	print $OUT join("\t", map {$data{$idx}{$_}} qw/ L C N /), "\n";
    }
    
} else {
    
    ## process Hskp_ncRNA bam
    open my $IN, '-|', "samtools view $hskpbam" or die "$0: Cannot 'samtools view' BAM file '$hskpbam': $!\n";
    print "TrimLen\tClass\tMulti\n";
    while (<$IN>) {
	my @fields = split /\t/, $_, 11;
	my $mapq = $fields[4];  # primary alignment score
	my ($sas) = ($fields[10] =~ /XS:i:(\d+)/);  # best secondary-alignment score
	my $multi = $sas >= $mapq-1 ? 1 : 0;  # if secondary score == primary score (or primary-1), assert that there were > 1 equally-good alignments
	my $class = $classes{ (split /\|/, $fields[2])[0] };  # type of ncRNA (or '*' if unaligned), as integer 0-8
	my $len = length($fields[9]);  # trimmed read length
	print $OUT "$len\t$class\t$multi\n";
    }
    close $IN;
}

close $OUT;
exit;
