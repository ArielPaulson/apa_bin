#!/usr/bin/perl

## this script automatically selects the most current GO database and runs FatiClone using it.
## accomplishes nothing important unless the database is brand new and unused
## if a new/unused database, FatiClone will cache the GO tree hash for future use

chomp(my @dbs = `/home/apa/local/bin/FatiClone --showdbs`);
foreach (@dbs) {
    if ($_ =~ /^(go_\d+)/) {
	$current = $1;
	last;
    }
}

print "Current db: $current\n";
mkdir '/home/apa/local/bin/GO_Tools_DBcache/tmp';

system "/home/apa/local/bin/FatiClone -f /home/apa/local/bin/GO_Tools_DBcache/FC_Dummy.txt -d $current -x 10090 -b genome -w /home/apa/local/bin/GO_Tools_DBcache/tmp";
exit;


