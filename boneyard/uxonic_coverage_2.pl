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
    my ($gene, $exons) = split /:/, $data[3];
    my $len = $data[2]-$data[1]+1;
    $exons =~ s/\|/,/g;      # convert | to , in uxon ID
    $length{$gene} += $len;   # uxonic length (1 BASED COORDS)
    $counts{$gene} = 0;  # ensure printable
    $covpct{$gene} = 0;  # ensure printable
    $allgenes{$gene}{$exons} = 1;  # store uxons by gene too
    $uxoncts{$gene}{$exons} = $uxoncvg{$gene}{$exons} = 0;  # ensure printable values
    $uxonlen{$gene}{$exons} = $len;
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
    $widths{$gene} += $width+1;    # coverageBed assumes 0-based
    $counts{$gene} += $count;
    $bpcov{$gene} += $nzb;
    $uxon =~ s/^$gene://;   # remove gene name from uxon ID
    $uxon =~ s/\|/,/g;      # convert | to , in uxon ID
    $uxoncts{$gene}{$uxon} = $count;
    $uxoncvg{$gene}{$uxon} = $nzb;
}
close IN;

$covpct{$_} = $bpcov{$_}/($widths{$_}||1) foreach keys %allgenes;

$outfile = "$bampath/gene_uxonic_coverage.txt" unless $outfile;
open OUT, "> $outfile" or die "Cannot write results to '$outfile': $!\n";
print OUT "Gene\tUxonLen\tCounts\tCovPct\tUxon|Len|Counts|CovBp\n";
foreach my $gene (sort keys %allgenes) {
    print OUT "$gene\t$length{$gene}\t$counts{$gene}\t$covpct{$gene}\t";
    print OUT "$_|$uxonlen{$gene}{$_}|$uxoncts{$gene}{$_}|$uxoncvg{$gene}{$_} " foreach sort keys %{ $allgenes{$gene} };
    print OUT "\n";
}
close OUT;



unlink $tempout;
chomp(my $now = `date`);
print "Complete: $now\n";
exit;

