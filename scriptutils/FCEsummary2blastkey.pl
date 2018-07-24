#!/usr/bin/perl

my @orfnames = qw/ +0 +1 +2 -0 -1 -2 /;
my %orfidx = ('+0',0, '+1',1, '+2',2, '-0',3, '-1',4, '-2',5);

open IN, $ARGV[0] or die "Cannot read '$ARGV[0]': $!\n";
open OUT, "> $ARGV[1]" or die "Cannot write '$ARGV[1]': $!\n";
while (<IN>) {
    next if $. < 3;
    $_ =~ s/[\n\r]+$//;
    my @data = split /\t/, $_;
    my ($ntlen, $potential) = @data[2,4];
    my $header = (split /\s/, $data[0])[0];
    my @orflens = @data[9,15,21,27,33,39];
    my @orfbest = (0,0,0,0,0,0);
    foreach my $best (split /,/, $data[1]) {
	$orfbest[$orfidx{$_}] = 1 foreach (split /,/, $data[1]);  # best orf names -> indexes -> @orfbest gets 1 not 0
    }
    print OUT "$header\t$ntlen\t$potential";
    print OUT "\t$orfnames[$_],$orflens[$_],$orfbest[$_]" foreach (0..5);
    print OUT "\n";
}
close IN;
close OUT;
exit;
