#!/usr/bin/perl

my $fasta = $ARGV[0];    # we are expecting the Trinity.fasta fasta
my ($fname) = ($fasta =~ /([^\/]+)$/);
my $outfile = "$fname.data.txt";

open IN, $fasta or die "Cannot read '$fasta': $!\n";
while (<IN>) {
    if ($_ =~ /^>/) {
	my ($contig, $len, $fpkm, $path) = ($_ =~ /^>(\S+) len=(\d+) ~FPKM=([\d.]+) path=\[(.*?)\]/);
	print "Failed to parse line: '$_'\n" unless $contig;
	push @data, "$contig\t$len\t$fpkm\t$path\n";
#	print "$contig | $len | $fpkm | $path\n";
    }
}
close IN;

#>comp50976_c1_seq2 len=111 ~FPKM=7 path=[0:0-19 20:20-29 30:30-53 141:54-110]

open OUT, "> $outfile" or die "Cannot write '$outfile': $!\n";
print OUT "Contig\tLength\tFPKM\tPath\n", @data;
close OUT;
exit;
