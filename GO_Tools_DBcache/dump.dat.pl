#!/usr/bin/env perl
use Data::Dumper;
use Storable qw/ retrieve /;

my %hash = %{ retrieve($ARGV[0]) };

open OUT,  "> $ARGV[0].dump.txt";
print OUT Dumper(\%hash), "\n";
close OUT;
exit;
