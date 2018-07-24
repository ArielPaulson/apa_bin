#!/usr/bin/env perl
require '/home/apa/apa_routines.pm';
use strict;


my ($infile, $conv, $outfile) = @ARGV;
my %conv;


my $CO = &open2('<', $conv);
while (<$CO>) {
    chomp;
    my ($old, $new) = split /\t/, $_;
    $conv{$old}{$new} = 1;
}
close $CO;


my $IN = &open2('<', $infile);
my $OUT = &open2('>', $outfile);
while (<$IN>) {
    chomp;
    if ($.==1) {
        print $OUT "$_\n";
    } else {
        my ($id, $rest) = split /\t/, $_, 2;
        print $OUT "$_\t$rest\n" foreach sort keys %{ $conv{$id} };
    }
}
close $IN;
close $OUT;
