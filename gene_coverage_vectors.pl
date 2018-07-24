#!/usr/bin/perl

my ($genebed, $bedgraph, $outfile) = @ARGV;
my ($gbname) = ($genebed =~ /([^\/]+)\.bed$/);
my ($bgname) = ($bedgraph =~ /([^\/]+)\.b[edgraph]{1,7}$/);
$outfile = "$bedgraph.$gbname.gene.vectors.txt" unless $outfile;
die "Bedgraph file '$bedgraph' does not exist!\n" unless -e $bedgraph;
die "Gene bed file '$genebed' does not exist!\n" unless -e $genebed;
my $tempfile = "temp.$$";

print "Intersecting...\n";
system "intersectBed -wao -a $genebed -b $bedgraph > $tempfile";

print "Reading...\n";
open IN, $tempfile or die "Couldn't read intersectBed output file '$tempfile': $!\n";
while (<IN>) {
    my @line = split /\t/, $_;
    my ($gstart, $gend, $gene, $strand, $start, $end, $height) = @line[1..3,5,7..9];
    my ($gene, $else) = split /:/, $gene;
    $data{$gene}{STR} = $strand;
    $data{$gene}{LEN} += $gend-$gstart+1;
    push @{ $data{$gene}{COORDS} }, ($gstart, $gend);
    if ($height ne '.') {
	$data{$gene}{MAP}{$_} += $height foreach ($start+1..$end);  # this may overrun gene boundaries
    }
}
close IN;
unlink $tempfile;

print "Writing $outfile...\n";
open OUT, "> $outfile" or die "Couldn't write output file '$outfile': $!\n";
foreach my $gene (sort keys %data) {
    my ($start, $end) = (sort {$a <=> $b} @{ $data{$gene}{COORDS} })[0,-1];
    my %temp = %{ $data{$gene}{MAP} };
    foreach (keys %temp) {
	delete $temp{$_} if ($_ < $start || $_ > $end);  # remove stuff beyond gene boundaries
    }
    my @vec = map { $temp{$_} } ($start..$end);   # ensure complete map
    @vec = reverse(@vec) if $data{$gene}{STR} eq '-';
    my $vec = join ',', @vec;
    print OUT "$gene\t$data{$gene}{LEN}\t$data{$gene}{STR}\t$vec\n";
}
close OUT;
exit;
