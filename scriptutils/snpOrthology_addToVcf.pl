#!/usr/bin/env perl
require '/home/apa/apa_routines.pm';
use Getopt::Long;
use Pod::Usage;
use strict;


## Adds snpOrthology data to VCF


## Inputs
my ($var, $invcf, $outvcf) = @ARGV;     # ___.variation_table.txt output of snpOrthology

my (%vdat, %alltrans, @cols, @VCF, @lines);

## Get Gene->Transcript Relations
my $VAR = &open2('R', $var, 'variation_table.txt file');
while (<$VAR>) {
    s/[\n\r]+$//;
    my @data = split /\t/, $_;
    if ($.==1) {
        @cols = @data;
    } else {
        my ($chr, $pos, $ref, $alt, $trans) = @data[0..4];
        my %tmp = map {($cols[$_]=>$data[$_])} (0..$#data);
        my $key = "$chr.$pos.$ref.$alt";
        $alltrans{$trans} = 1;
        $vdat{$key}{T}{$trans} = \%tmp;
    }
}
close $VAR;


## Read VCF
my $IN = &open2('R', $invcf, 'input vcf');
while (<$IN>) {
    s/[\n\r]+$//;
    
    if (/^#/) {
        
        ## Headers
        push @VCF, "$_\n";
        
    } else {
        
        ## Records
        $lines[0]++;
        my @fields = split /\t/, $_;
        my ($chr, $pos, $id, $refs, $alts) = @fields[0..4];
        my %gt;
        
        ## Look for snpEff data in INFO field
        foreach my $item (split /;/, $fields[7]) {
            my @keyval = split /=/, $item;
            if ($keyval[0] eq 'EFF') {
                ## Old snpEff effects
                die "$0: gene-transcript relations are required for older EFF-type snpEff annotations!\n"; # if $gene_no_gtr;
                foreach my $entry (split /,/, $keyval[1]) {
                    my $tid = (split /\|/, $entry)[9];
                    $gt{$tid} = 1;
                }
            } elsif ($keyval[0] eq 'ANN') {
                ## New snpEff effects
                foreach my $entry (split /,/, $keyval[1]) {
                    my ($eff, $gid, $tid) = (split /\|/, $entry)[1,4,6];
                    $gt{$tid} = 1 if $eff =~ /missense/i;
                }
            }
        }
        
        ## Get unique genes/transcripts referred to
        if (%gt) {
            my @altortho;
            $lines[1]++;
            foreach my $ref (split /,/, $refs) {
                foreach my $alt (split /,/, $alts) {
                    my $key = join('.', $chr, $pos, $ref, $alt);
                    if (exists $vdat{$key}) {
                        $vdat{$key}{N}++;
                        foreach my $trans (keys %gt) {
                            my $trans2 = $trans;
                            $trans2 =~ s/\.\d+$// unless exists $vdat{$key}{T}{$trans2};
                            next unless exists $vdat{$key}{T}{$trans2};
                            my ($Talt, $AltObs, $RefClass, $AltClass, $OrthoR) = map { $vdat{$key}{T}{$trans2}{$_} } qw/ AltNT.T Alt.Obs Ref.Class Alt.Class OrthoR /;
                            push @altortho, "$alt|$Talt|$trans2|$AltObs|$RefClass|".($RefClass-$AltClass).".$OrthoR";
                        }
                    }
                }     
            }
            if (@altortho) {
                $lines[2]++;
                $fields[7] .= ";SnpOrtho=".join(',', @altortho);
            }
        }
        
        
        push @VCF, join("\t", @fields)."\n";
        
    }
}
close $IN;


## Write VCF
my $OUT = &open2('W', $outvcf, 'output VCF');
print $OUT @VCF;
close $OUT;

my @matched;
my $Nkeys = scalar keys %vdat;
foreach my $key (keys %vdat) {
    $vdat{$key}{N} ? $matched[1]++ : $matched[0]++;
}

## Reporting and exit
print "$lines[0] records total.\n$lines[1] records with missense events.\n$lines[2] records had events added to the INFO field.\n";
print "$Nkeys \$var rows.\n$matched[1] rows matched in VCF.\n$matched[0] rows unmatched.\n";
print "$0 $invcf complete!\n";
exit;


