#!/usr/bin/perl

## Takes a BAM file and an uxon BED file and returns (approximate) genewise read counts. 
## 20110818 by apa.  

## v2 requires no output file, instead writes generic-named file to $bam directory

my ($uxonbed, $bam, $outfile) = @ARGV;   # $uxonbed e.g. from my bowtie index sets
my ($bampath) = ($bam =~ /^(.*)\/[^\/]+$/);
$bampath = '.' unless $bampath;   # bam file might not have explicit path
my $tempout = "uxon_covg_temp_$$.txt";

chomp(my $now = `date`);
print "Reading $uxonbed : $now\n";
open IN, $uxonbed or die "Cannot open uxon bed file '$uxonbed': $!\n";
while (<IN>) {
    next if $_ =~ /^track/;
    $_ =~ s/[\n\r]+$//;
    my @data = split /\t/, $_;
    my ($gene, $else) = split /:/, $data[3];
    $length{$gene} += $data[2]-$data[1]+1;   # uxonic length (1 BASED COORDS)
    $allgenes{$gene} = 1;
    $counts{$gene} = 0;  # ensure printable
    $covpct{$gene} = 0;  # ensure printable
}
close IN;

my $com = "coverageBed -split -abam $bam -b $uxonbed > $tempout";
chomp(my $now = `date`);
print "$com : $now\n";
system $com;

chomp(my $now = `date`);
print "Post-processing $tempout : $now\n";
open IN, $tempout or die "Cannot open uxon coverage file '$tempout': $!\n";
while (<IN>) {
    $_ =~ s/[\n\r]+$//;
    my @data = split /\t/, $_;
    my ($uxon, $count, $nzb, $width) = @data[3,6..8];
    my ($gene, $else) = split /:/, $uxon;
    $widths{$gene} += $width + 1;    # coverageBed assumes 0-based
    $counts{$gene} += $count;
    $bpcov{$gene} += $nzb;
}
close IN;

$covpct{$_} = $bpcov{$_}/($widths{$_}||1) foreach keys %allgenes;

$outfile = "$bampath/gene_uxonic_coverage.txt" unless $outfile;
open OUT, "> $outfile" or die "Cannot write results to '$outfile': $!\n";
print OUT "Gene\tUxonLen\tCounts\tCovPct\n";
print OUT "$_\t$length{$_}\t$counts{$_}\t$covpct{$_}\n" foreach sort keys %allgenes;
close OUT;

unlink $tempout;
chomp(my $now = `date`);
print "Complete: $now\n";
exit;

