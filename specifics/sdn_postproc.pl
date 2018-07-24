#!/usr/bin/perl
use Cwd;

my $wd = cwd();

foreach my $N (180,300) {
    foreach my $k (25,50,75) {
	chdir "$wd/inornata_$N/$k";
	my $tag = "inornata_$N.$k";
	system "grep '>' $tag.contig > $tag.contig.headers";
	system "cut -d' ' -f3 $tag.contig.headers > $tag.contig.lengths";
    }
}
print "Complete!\n";
exit;
