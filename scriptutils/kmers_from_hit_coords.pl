#!/usr/bin/perl
use Getopt::Long;

# this script was developed to take a bed file of probe blat-hit coords
# and return a bed file of coords of all kmers which incorporated the blat hits.
# This to simulate what reads of length k we might expect from doing region-capture.

my ($bed, $k, $overlap, $zero);  # zero to indicate 0-based coords (NOT READY)
GetOptions("f=s" => \$bed, "k=i" => \$k, "o=i" => \$overlap, "zero" => \$zero);

my $overhang = $k - $overlap;
die "Target overlap '$overlap' cannot exceed kmer size '$k'!\n" if $overhang < 0;

open IN, $bed or die "Couldn't open file '$bed': $!\n";
open OUT, "> ${k}mers.$bed";
while (<IN>) {
    next if $_ =~ /^track/ || $_ =~ /^#/;
    $_ =~ s/[\n\r]+$//;
    my @data = split /\t/, $_;
    my ($chr, $start, $end, $id) = @data[0..3];
    my $len = $end-$start+1;
    my $kstart = $start-$overhang;
    my $kend = $end+$overhang;
#    if ($k > $len) {
#	($kstart, $kend) = ($start-($k-$len), $start);
#    } elsif ($k < $len) {
#	($kstart, $kend) = ($start, $length-$k+1);
#    } elsif ($k == $len) {
#	($kstart, $kend) = ($start, $end);
#    }
    my $n;
    foreach my $i ($kstart..$kend) {
	$n++;
	print OUT "$chr\t$i\t",($i+$k-1),"\t$id:$n\n";
    }
}
close IN;
close OUT;
exit;
