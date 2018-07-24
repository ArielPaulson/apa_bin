#!/usr/bin/perl

my ($bam) = @ARGV;
my ($name) = ($bam =~ /([^\/]+)\.bam$/i);
my $outfile = "$name.read.uniqueness.txt";
die "BAM file '$bam' does not exist!\n" unless -e $bam;
my $temp = "temp.$$";

system "samtools view $bam | cut -f1,10 > $temp";

open IN, $temp;
while (<IN>) {
    chomp;
    my ($read, $seq) = split /\t/, $_, 2;
    $aligns++;
    $reads{$read}++;
    $seqs{$seq}++;
}
close IN;
unlink $temp;

$Nreads = scalar keys %reads;
$Nseqs = scalar keys %seqs;

open OUT, "> $outfile" or die "Couldn't write output to '$outfile': $!\n";
print OUT "ALIGNS\t$aligns\n";
print OUT "READS\t$Nreads\n";
print OUT "SEQUENCES\t$Nseqs\n";
close OUT;
exit;

