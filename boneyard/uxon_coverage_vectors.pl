#!/usr/bin/perl

my ($uxonbed, $bedgraph, $outfile) = @ARGV;
my ($ubname) = ($uxonbed =~ /([^\/]+)\.bed$/);
my ($bgname) = ($bedgraph =~ /([^\/]+)\.b[edgraph]{1,7}$/);
$outfile = "$bedgraph.$ubname.uxon.vectors.txt" unless $outfile;
die "Bedgraph file '$bedgraph' does not exist!\n" unless -e $bedgraph;
die "Uxon bed file '$uxonbed' does not exist!\n" unless -e $uxonbed;
my $tempfile = "temp.$$";

print "Intersecting...\n";
system "intersectBed -wao -a $uxonbed -b $bedgraph > $tempfile";

print "Reading...\n";
open IN, $tempfile or die "Couldn't read intersectBed output file '$tempfile': $!\n";
while (<IN>) {
    my @line = split /\t/, $_;
    my ($ustart, $uend, $uxon, $strand, $start, $end, $height) = @line[1..3,5,7..9];
    my ($gene, $else) = split /:/, $uxon;
    $data{$gene}{STR} = $strand;
    $data{$gene}{LEN} += $uend-$ustart+1;
    $data{$gene}{COORDS}{"$ustart\t$uend"} = 1;
    if ($height ne '.') {
	$data{$gene}{MAP}{$_} += $height foreach ($start+1..$end);  # this may overrun gene boundaries
    }
}
close IN;
unlink $tempfile;

print "Writing $outfile...\n";
open OUT, "> $outfile" or die "Couldn't write output file '$outfile': $!\n";
foreach my $gene (sort keys %data) {
    my @vec;
    foreach my $uxon (sort {$a <=> $b} keys %{ $data{$gene}{COORDS} }) {
	my ($ustart, $uend) = split /\t/, $uxon;
	push @vec, $data{$gene}{MAP}{$_} foreach ($ustart..$uend);   # ensure complete map
    }
    @vec = reverse(@vec) if $data{$gene}{STR} eq '-';
    my $vec = join ',', @vec;
    print OUT "$gene\t$data{$gene}{LEN}\t$data{$gene}{STR}\t$vec\n";
}
close OUT;
exit;
