#!/usr/bin/perl

my @letters = (A..Z);
open IN, $ARGV[0] or die "Couldn't read file '$ARGV[0]': $!\n";
open OUT, "> $ARGV[0].bed" or die "Couldn't write file '$ARGV[0].bed': $!\n";
while (<IN>) {
    next unless $_ =~ /^\d/;  # no headers
    $_ =~ s/[\n\r]+$//;
    my @data = split /\t/, $_;
    my ($score, $strand, $id, $chr, $start, $end, $lengths, $blocks) = @data[0,8,9,13,15,16,18,20];
    $lengths =~ s/,$//;
    $blocks =~ s/,$//;
    my @blocks = (split /,/, $blocks);
    if (scalar @blocks == 1) {  # gapless
	print OUT "$chr\t$start\t$end\t$id\t$score\t$strand\n";
	$gapless++;
    } else {
	my @lengths = (split /,/, $lengths);
	foreach my $i (0..$#blocks) {
	    my ($bstart, $bend, $letter) = ($blocks[$i], $blocks[$i]+$lengths[$i], $letters[$i]);
	    print OUT "$chr\t$bstart\t$bend\t$id:$letter\t$score\t$strand\n";
	}
	$gapped++;
    }
}
close IN;
close OUT;
print "$ARGV[0]: $gapless gapless | gapped $gapped\n";
exit;
