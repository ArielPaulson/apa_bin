#!/usr/bin/env perl
use Net::FTP;
use strict;

#die;	# ALREADY RUN

#chdir '/n/data1/genomes/Sanger_Mouse_Genomes/' or die;

my $ftp = Net::FTP->new('ftp.sanger.ac.uk') or die "Cannot connect: $!\n";
$ftp->login('anonymous','apa@stowers.org');

$ftp->binary();

$ftp->cwd('/pub/mouse_genomes/current_consensus/');
$ftp->get('README.txt');
my @dirs = $ftp->ls('*_Mouse_Genome');
foreach my $dir (@dirs) {
	mkdir $dir unless -d $dir;
	my @files = $ftp->ls($dir);
	foreach my $file (@files) {
		print "Getting $dir/$file...\n";
		$ftp->get("$dir/$file", "$dir/$file");
	}
}
$ftp->quit or die "Couldn't close connection cleanly: $!\n";	
exit;
