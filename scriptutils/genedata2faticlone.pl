#!/usr/bin/env perl

my ($org, $date) = @ARGV;  # e.g. dme 201505
my $prefix = "$org/$date/${org}_$date";
my $input = "${prefix}_gene_data.txt";
my $output = "${prefix}_FatiClone_DB.txt";

my %data;

open my $IN, '<', $input or die "$0: cannot read input file '$input: $!\n";
open my $OUT, '>', $output;
my $id_col = 2;  # third column (starting from 0) -- this seems to hold across orgs
my $path_col;  # this does not
while (<$IN>) {
    chomp;
    if ($. == 1) {
	my @colnames = split /\t/, $_;
	foreach my $i (0..$#colnames) {
	    if ($colnames[$i] eq 'Pathways') {
		$path_col = $i;
		last;
	    }
	}
    }
    my ($id, $paths) = (split /\t/, $_)[$id_col,$path_col];
    next unless $paths;
    foreach my $entry (split /; /, $paths) {
	my ($pathid, $pathname) = split / /, $entry, 2;
	print $OUT "$id\tKEGG\t$pathid\t$pathname\n";
    }
}
close $IN;
close $OUT;
exit;



