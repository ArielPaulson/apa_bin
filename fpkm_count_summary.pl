#!/usr/bin/env perl

my $gbed = '/n/data1/genomes/bowtie-index/hg19/Ens_66/hg19.Ens_66.genes.bed';
my $tbed = '/n/data1/genomes/bowtie-index/hg19/Ens_66/hg19.Ens_66.trans.bed';
my $nbed = '/n/projects/apa/stuff/Noncode/human/human.bed';

print "Reading gene definitions...\n";
open IN, $gbed or die;
while (<IN>) {
    next if $. == 1;
    my @data = split /\t/, $_;
    $allids{$_}{$data[3]} = 1 foreach qw/ genes counts /;
}
print "$. genes\n";
close IN;

open IN, $tbed or die;
while (<IN>) {
    next if $. == 1;
    my @data = split /\t/, $_;
    $allids{isoforms}{$data[3]} = 1;
}
print "$. isoforms\n";
close IN;

open IN, $nbed or die;
while (<IN>) {
    my @data = split /\t/, $_;
    $allids{$_}{$data[3]} = 1 foreach qw/ genes counts /;
    $allids{isoforms}{"$data[3].1"} = 1;
}
print "$. noncoding\n";
close IN;

print "Reading FPKMs...\n";
my ($FPKM, $COVG);
#chdir 'data';
foreach my $sample (glob "*") {
    next unless -d $sample;
    next unless -e "$sample/$sample.bam";  # otherwise not a sample directory
    chdir $sample;
    print "$sample\n";
    foreach my $level (qw/ genes isoforms /) {
	open IN, "$level.fpkm_tracking" or print "Couldn't find '$level.fpkm_tracking': $!\n";
	while (<IN>) {
	    $_ =~ s/[\n\r]+$//;
	    my @data = split /\t/, $_;
	    if ($. == 1) {
		foreach my $i (0..$#data) {
		    $FPKM = $i if $data[$i] eq 'FPKM';
		    $COVG = $i if $data[$i] eq 'coverage';
		}
		die "FPKM column not identified!  Header: '$_'\n" unless $FPKM;
		die "Coverage column not identified!  Header: '$_'\n" unless $COVG;
	    } else {
		my ($ID, $covg, $fpkm) = @data[0,$COVG,$FPKM];
		$master{$sample}{$level}{$ID}{FPKM} += $fpkm;
		$master{$sample}{$level}{$ID}{COVG} += $covg;
		$master{$sample}{$level}{$ID}{N}++;
	    }
	}
	close IN;
    }
    if (-e "gene_uxonic_coverage.txt") {
	open IN, "gene_uxonic_coverage.txt" or print "Couldn't find 'gene_uxonic_coverage.txt': $!\n";
	while (<IN>) {
	    $_ =~ s/[\n\r]+$//;
	    my ($ID, $uxlen, $count, $covpct) = split /\t/, $_;
	    $master{$sample}{counts}{$ID}{FPKM} += $count;
	    $master{$sample}{counts}{$ID}{COVG} += $covpct;
	    $master{$sample}{counts}{$ID}{N}++;
	    $uxlens{$ID} = $uxlen;
	}
	close IN;
    }
    chdir '..';
}
#chdir '..';

my $tag = 'all';
#my $tag = 'u50';
print "Writing...\n";
my %labels = ('genes' => 'Gene', 'isoforms' => 'Transcript', 'counts' => 'Gene');
foreach my $level (qw/ genes isoforms counts /) {
    if ($level eq 'counts') {
	open OUT, "> $tag.genes.$level.txt";
    } else {
	open OUT, "> $tag.$level.fpkms.txt";
    }
    print OUT $labels{$level};
    print OUT "\tUxonLen" if $level eq 'counts';
    foreach my $sample (sort keys %master) {
	print OUT "\t$sample";
	print OUT "\tCOVG" if $level ne 'genes';
    }
    print OUT "\n";
    foreach my $ID (sort keys %{ $allids{$level} }) {
	print OUT $ID;
	print OUT "\t$uxlens{$ID}" if $level eq 'counts';
	foreach my $sample (sort keys %master) {
	    print OUT "\t$master{$sample}{$level}{$ID}{FPKM}";
	    if ($level ne 'genes') {
		my $covg = $master{$sample}{$level}{$ID}{N} ? sprintf("%0.2f", 100*$master{$sample}{$level}{$ID}{COVG}/$master{$sample}{$level}{$ID}{N}) : 0;
		print OUT "\t$covg";
	    }
	}
	print OUT "\n";
    }
    close OUT;
}
print "Complete!\n";
exit;
