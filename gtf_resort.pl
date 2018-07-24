#!/usr/bin/perl
use strict;

my ($in, $out) = @ARGV;

my (@header, @chrs, %data, %gt, %already, $i);

open IN, '<', $in or die;
while (<IN>) {
    if (/^#/) {
        push @header, $_;
        $i++;
    } else {
        my @x = split /\t/, $_;
        my ($gene, $trans);
        ($gene) = ($x[8] =~ /gene_id "([^"]+)"/);
        ($trans) = ($x[8] =~ /transcript_id "([^"]+)"/);
        die "$0: line $.: no gene id!\n$_\n" unless $gene;
        die "$0: line $.: no trans id!\n$_\n" if !$trans && $x[2] ne 'gene';
        if ($x[2] eq 'gene') {
            $gt{G}{$gene} = $_;
            $i++;
        } elsif ($x[2] eq 'transcript') {
            $gt{T}{$trans} = $_;
            $i++;
        } else {
            push @chrs, $x[0] unless exists $data{$x[0]};
            push @{ $data{$x[0]}{$x[3]}{$x[4]}{$x[2]} }, [$_, $gene, $trans];
            $i++;
        }
    }
}
close IN;

open OUT, '>', $out or die;
print OUT @header;
my $o;
foreach my $chr (@chrs) {
    foreach my $start (sort {$a <=> $b} keys %{ $data{$chr} }) {   # ASCENDING
        foreach my $end (sort {$b <=> $a} keys %{ $data{$chr}{$start} }) {   # DESCENDING
            #my $nfeat = scalar keys %{ $data{$chr}{$start}{$end} };
            foreach my $feat (sort {$b cmp $a} keys %{ $data{$chr}{$start}{$end} }) {   # DESCENDING
                #print "$start-$end: $feat\n" if $nfeat > 1;
                foreach my $entry (@{ $data{$chr}{$start}{$end}{$feat} }) {
                    my ($line, $gene, $trans) = @$entry;
                    warn "Missing G or T!\n$line\n" unless $gene && $trans;
                    print OUT $gt{G}{$gene} if exists $gt{G}{$gene} && ! exists $already{G}{$gene};
                    print OUT $gt{T}{$trans} if exists $gt{T}{$trans} && ! exists $already{T}{$trans};
                    $already{G}{$gene} = 1;
                    $already{T}{$trans} = 1;
                    print OUT $line;
                    #print "$chr, $start, $end, $feat:  $line";
                    $o++;
                    #die if $o > 20;
                }
            }
        }
    }
}
close OUT;
print "$i, $o\n";
exit;
