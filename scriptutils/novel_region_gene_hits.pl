#!/usr/bin/perl

die unless scalar @ARGV >= 4;  # must specify at least 4 files
my $out = $ARGV[4];
$out = 'novel_region_gene_hits.txt' unless $out;

my %inputs = (
    'exons' => { 
        72 => $ARGV[0],  # 'novel.probes.72.8mm.psl.bed_mm9.Ens_63.exons.bed_N_N.txt',
        100 => $ARGV[1],  # 'novel.probes.100.8mm.psl.bed_mm9.Ens_63.exons.bed_N_N.txt'
    },
    'genes' => { 
        72 => $ARGV[2],  # 'novel.probes.72.8mm.psl.bed_mm9.Ens_63.genes.bed_N_N.txt',
        100 => $ARGV[3],  # 'novel.probes.100.8mm.psl.bed_mm9.Ens_63.genes.bed_N_N.txt'
    }
);

push @output, "K\tSet\tIntergenic\tGenic.Exonic\tGenes.Exons\n";

foreach my $set (qw/ exons genes /) {
    foreach my $size (72,100) {
        my ($m, $u, $g);
        open IN, "match_stats_$inputs{$set}{$size}" or die;
        while (<IN>) {
            $_ =~ s/[\n\r]+$//;
            if ($_ =~ /^Elements matched/) {
                my @data = split /\s+/, $_;
                $m = $data[2];
                $g = $data[3];
            } elsif ($_ =~ /^Elements unmatched/) {
                my @data = split /\s+/, $_;
                $u = $data[2];
            }
        }
        close IN;
        push @output, "$size\t$set\t$u\t$m\t$g\n";
    }
}

push @output, "\nK\tSet\tGene\tBed\n";

foreach my $set (qw/ exons genes /) {
    foreach my $size (72,100) {
        open IN, "matched_$inputs{$set}{$size}" or die;
        while (<IN>) {
            next if $. == 1;
            $_ =~ s/[\n\r]+$//;
            my @data = split /\s+/, $_;
            my ($chr, $start, $end, $exon, $gene) = @data[0..3,7];
            $exon =~ s/^\d+://;
            my $coord = join ':', ($chr, $start, $end, $exon);
            push @output, "$size\t$set\t$gene\t$coord\n";
        }
        close IN;
    }
}

open OUT, "> $out";
print OUT @output;
close OUT;
exit;


