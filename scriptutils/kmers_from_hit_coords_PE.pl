#!/usr/bin/perl
use Getopt::Long;

# this script was developed to take a bed file of probe blat-hit coords
# and return a bed file of coords of all kmers which incorporated the blat hits.
# This to simulate what reads of length k we might expect from doing region-capture.

my ($bed, $k, $overlap, $insert, $zero);  # zero to indicate 0-based coords (NOT READY)
GetOptions("f=s" => \$bed, "k=i" => \$k, "o=i" => \$overlap, "i=i" => \$insert, "zero" => \$zero);

my $overhang = $k - $overlap;
die "Target overlap '$overlap' cannot exceed kmer size '$k'!\n" if $overhang < 0;

open IN, $bed or die "Couldn't open file '$bed': $!\n";
open OUT1, "> ${k}mers.PE1_${insert}_$bed";
open OUT2, "> ${k}mers.PE2_${insert}_$bed";
while (<IN>) {
    next if $_ =~ /^track/ || $_ =~ /^#/;
    $_ =~ s/[\n\r]+$//;
    my @data = split /\t/, $_;
    my ($chr, $start, $end, $id) = @data[0..3];
    my $len = $end-$start+1;
    my $kstart1 = $start-$overhang-$insert-$k;
    my $kend1 = $end-$overlap-$insert-$k;  # last End2 start = $end-$overlap, thus last End1 start = that-$insert-$k
    my $kstart2 = $start-$overhang;
    my $kend2 = $end-$overlap;
    my $n;
    foreach my $i ($kstart1..$kend1) {
	my $j = $i + $k + $insert;
	$n++;
	print OUT1 "$chr\t$i\t",($i+$k-1),"\t$id:$n:1\n";
	print OUT2 "$chr\t$j\t",($j+$k-1),"\t$id:$n:2\n";
    }
    foreach my $i ($kstart2..$kend2) {
	my $j = $i + $k + $insert;
	$n++;
	print OUT1 "$chr\t$i\t",($i+$k-1),"\t$id:$n:1\n";
	print OUT2 "$chr\t$j\t",($j+$k-1),"\t$id:$n:2\n";
    }
}
close IN;
close OUT1;
close OUT2;
exit;
