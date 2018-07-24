#!/usr/bin/perl
use LWP::UserAgent;
use Getopt::Long;

my ($orgcode, $outdir);
GetOptions("c=s" => \$orgcode, "o=s" => \$outdir);

$orgcode = "\L$orgcode";
my $ua = LWP::UserAgent->new;

#system "rm -Rf $outdir";
system "mkdir -p $outdir";
die "Cannot create output directory '$outdir' in this location!\n" unless -d $outdir;
chdir $outdir;

chomp(my $begin = `date`);
print "Beginning: $begin\n";

my ($page1, $page2, $dump);

## Build Gene Queries

$dump = "gene_list_dbget.txt";
if (-s $dump) {
	print "Reading gene list query dump...\n";
	open IN, $dump;
	chomp(my @page2 = (<IN>));
	$page2 = join '', @page2;
	close IN;
} else {
	print "Querying gene list...\n";
	my $url = 'http://www.genome.jp/dbget-bin/www_bfind_sub?max_hit=100000&dbkey=genes&keywords='.$orgcode.'%3A&mode=bfind';
	print " GET: $url\n";
	my $request = HTTP::Request->new(GET => $url);
	my $response = $ua->request($request);
	if ($response->is_success) {
		&query_dump($response->content, $dump);
	} else {
		die "dbget 1 failed: ", $response->status_line, "\n";
	}
	sleep 1;
}

while ($page2 =~ m!<a href="/dbget-bin/www_bget\?($orgcode:\d+)">\1</a><br><div style="margin-left:2em">(.*?)</div>!mg) {
#	print "$1 = $2\n";
	$genes{$1} = $2;
	push @list, $1;
	$hits++;
	if ($hits == 200) {
		$list = join '+', @list;
		push @geneQueries, "http://www.genome.jp/dbget-bin/www_bget?$list";
		@list = ();
		$hits = 0;
	}
}
$list = join '+', @list;
push @geneQueries, "http://www.genome.jp/dbget-bin/www_bget?$list";


## Run/Parse Gene Queries

my $N = $#geneQueries+1;
foreach my $i (1..$N) {
	$dump = "gene_data_dbget_$i.txt";
	if (-s $dump) {
		print "Reading gene data query dump $i...\n";
		open IN, $dump;
		chomp(my @page2 = (<IN>));
		$page2 = join '', @page2;
		close IN;
	} else {
		print "Querying gene data $i/$N...\n";
		print " GET: ",$geneQueries[$i-1],"\n";
		my $request = HTTP::Request->new(GET => $geneQueries[$i-1]);
		my $response = $ua->request($request);
		if ($response->is_success) {
			&query_dump($response->content, $dump);
		} else {
			print "dbget 2 (i=$i) failed: ", $response->status_line, "\n";
		}
		sleep 1;
	}
	
	foreach my $line (split /\s/, $page2) {
		if ($line =~ /^ENTRY\s+(\S+)/) {
			$id = $1;
			$ids{$id} = 1;
		} elsif ($line =~ /^NAME\s+(\S+)/) {
			$name{$id} = $1;
		} elsif ($line =~ /^DEFINITION\s+(.*)/) {
			$def{$id} = $1;
		} elsif ($line =~ /^ORTHOLOGY\s+(.*)/) {
			push @{ $ortho{$id} }, $1;
			$hash = 'ortho';
		} elsif ($line =~ /^PATHWAY\s+(.*)/) {
			push @{ $path{$id} }, $1;
			$hash = 'path';
			$allpaths{$1} = 1;
		} elsif ($line =~ /^CLASS\s+(.*)/) {
			push @{ $class{$id} }, $1;
			$hash = 'class';
		} elsif ($line =~ /^MOTIF\s+(.*)/) {
			push @{ $motif{$id} }, $1;
			$hash = 'motif';
		} elsif ($line =~ /^DBLINKS\s+(.*)/) {
			push @{ $db{$id} }, $1;
			$hash = 'db';
		} elsif ($line =~ /^CODON/) {
			$hash = 0;
		} elsif ($line =~ /^AASEQ/) {
			$hash = 0;
		} elsif ($line =~ /^NTSEQ\s+(.*)/) {
			$seqlen{$id} = $1;
			@{ $seq{$id} } = ();
			$hash = 'seq';
		} elsif ($line =~ /^\/\/\//) {
			$hash = 0;
		} else {
			$line =~ s/^\s+//;
			push @{ $$hash{$id} }, $line if ($hash && $line =~ /\S/);
		}
	}
}

	
## Write Genewise Datasets

open KEY, "> key_$file";
print KEY "ID\tNAME\tDEFINITION\tORTHOLOGY\tPATHWAY\tCLASS\tMOTIF\tDBLINKS\n";
open FASTA, "> genes.fa";
open GENE, "> id_gene.txt";
open ORTHO, "> id_ortho.txt";
open PATH, "> id_path.txt";
open CLASS, "> id_class.txt";
open MOTIF, "> id_motif.txt";
open DBLS, "> id_dblinks.txt";
foreach $id (keys %ids) {
	my $Fsequence = join "\n", @{ $seq{$id} };
	my $Dsequence = join '', @{ $seq{$id} };
	$Fsequence =~ s/(.*)/\U$1/g;
	$Dsequence =~ s/(.*)/\U$1/g;
	print FASTA ">$id|$seqlen{$id}\n$Fsequence\n";
	my ($flag, $cumulative);
	foreach $line (@{ $class{$id} }) {					# this loop completely relies on the array being in input-order (FIFO)
		if ($flag == 0 && $line =~ /\[PATH:/) {			# complete entry
			push @{ $reclass{$id} }, $line;
		} elsif ($flag == 1 && $line =~ /\[PATH:/) {	# reached end-of-entry fragment (all entries must end with a "[PATH:{blah}]" statement)
			push @{ $reclass{$id} }, $cumulative.$line;
			($flag, $cumulative) = (0, '');
		} else {				# no "[PATH:" 
			$flag = 1;			# reached end of line, but not entry.
			$cumulative .= "$line ";
		}
	}
	my ($flag, $cumulative, $block);
	foreach $line (@{ $ortho{$id} }) {				# this loop completely relies on the array being in input-order (FIFO)
		if ($flag == 0 && $line =~ /\[EC:/) {		# complete entry
			push @{ $reortho{$id} }, $line;
			$block = 1;
		} elsif ($flag == 1 && $line =~ /^KO:/) {	# no ending "[EC:{blah}]"; reached start of new entry
			push @{ $reortho{$id} }, $cumulative;
			($flag, $cumulative) = (0, '');
			if ($line =~ /\[EC:/) {		# now test the new line: obviously complete entry or no?
				push @{ $reortho{$id} }, $line;
			} else {
				$flag = 1;				# reached end of line, but maybe not entry.
				$cumulative .= "$line ";
			}
		} elsif ($flag == 1 && $line =~ /\[EC:/) {	# reached end-of-entry fragment ("[EC:{blah}]", if present, denotes end)
			push @{ $reortho{$id} }, $cumulative.$line;
			($flag, $cumulative) = (0, '');
		} else {				# no "[EC:", regardless of flag 
			$flag = 1;			# reached end of line, but maybe not entry.
			$cumulative .= "$line ";
		}
	}
	push @{ $reortho{$id} }, $cumulative if ($flag == 0 && scalar @{ $ortho{$id} } == 1 && $block != 1);
	push @{ $reortho{$id} }, $cumulative if ($flag == 1 && $cumulative ne '');
	my ($flag, $mstring, $pcount);
	foreach $line (@{ $motif{$id} }) {	# this loop completely relies on the array being in input-order (FIFO)
		if ($line =~ /(Pfam|PROSITE): /) {
			push @{ $remotif{$id} }, $mstring if (defined $mstring);	# so, not on first entry
			$pcount++;
			($flag, $mstring) = (0, $line);
		} else {
			$flag = 1;
			$mstring .= " $line";
		}
	}
	push @{ $remotif{$id} }, $mstring if ($flag == 1 || scalar @{ $motif{$id} } == 1);
	push @{ $remotif{$id} }, $mstring if ($flag == 0 && $pcount > 1);
	my $orthos = join ' // ', @{ $reortho{$id} };
	my $paths = join ' // ', @{ $path{$id} };
	my $classes = join ' // ', @{ $reclass{$id} };
	my $motifs = join ' // ', @{ $remotif{$id} };
	my $dbs = join ' // ', @{ $db{$id} };
	print KEY "$id\t$name{$id}\t$def{$id}\t$orthos\t$paths\t$classes\t$motifs\t$dbs\n";
	print GENE "$id\t$name{$id}\t$def{$id}\t$Dsequence\n";
	print ORTHO "$id\t$_\n" foreach @{ $reortho{$id} };
	print PATH "$id\t$_\n" foreach @{ $path{$id} };
	print CLASS "$id\t$_\n" foreach @{ $reclass{$id} };
	print MOTIF "$id\t$_\n" foreach @{ $remotif{$id} };
	print DBLS "$id\t$_\n" foreach @{ $db{$id} };
}
close $_ foreach qw/ KEY GENE ORTHO PATH CLASS MOTIF DBLS FASTA /;


## Pathway Queries

my (%lists, @queries, @output);
my ($listcount, $capture, $pathway, $fieldno);
my $listno = 1;

my (@pathQueries, @list, $hits);
foreach my $path (keys %allpaths) {
	push @list, $path;
	$hits++;
	if ($hits == 200) {
		$list = join '+', @list;
		push @pathQueries, "http://www.genome.jp/dbget-bin/www_bget?pathway:$list";
		@list = ();
		$hits = 0;
	}
}
$list = join '+', @list;
push @pathQueries, "http://www.genome.jp/dbget-bin/www_bget?pathway:$list";


my ($pathway, %targets);
my $N = $#pathQueries+1;
foreach my $i (1..$N) {
	$dump = "pathway_data_dbget_$i.txt";
	if (-s $dump) {
		print "Reading pathway data query dump $i...\n";
		open IN, $dump;
		$page1 = join '', (<IN>);
		close IN;
	} else {
		print "Querying pathway data $i/$N:...\n";
		print " GET: ",$pathQueries[$i-1],"\n";
		my $request = HTTP::Request->new( GET => $pathQueries[$i-1] );
		my $response = $ua->request($request);
		if ($response->is_success) {
			&query_dump($response->content, $dump);
		} else {
			print "dbget 3 (i=$i) failed: ", $response->status_line, "\n";
		}
		sleep 1;
	}
	
	foreach my $line (split /\n/, $page1) {
		if ($line =~ /<nobr>(.*)<\/nobr><\/th>/) {
			$capture = $1;
		} elsif ($capture eq 'Entry') {
			($pathway) = ($line =~ /<code>(.*)/);
			$pathway =~ s/&nbsp;//g;
			$pathway =~ s/Pathway//g;
			$pathway =~ s/<br>//g;
			$capture = undef;
	#        print "Pathway: $pathway\n";
		} elsif ($capture eq 'Name') {
			my ($prename) = ($line =~ />([^<>]+)<br>/);
			$targets{$pathway}{Name} = (split / - /, $prename)[0];
			$capture = undef;
	#        print "  Name: $targets{$pathway}{Name}\n";
		} elsif ($capture eq 'Description') {
			$line =~ s/\[MD:[^\[\]]+\]//g;      # kill these things first
			if ($line =~ /^<.*>([^<>]+)<br>$/) {     # first/last line
				$targets{$pathway}{Description} .= $1;
			} elsif ($line =~ /^<.*>([^<>]+)$/) {    # first line
				$targets{$pathway}{Description} .= $1;
			} elsif ($line =~ /^([^<>]+)<br>$/) {    # last line
				$targets{$pathway}{Description} .= $1;
			} elsif ($line !~ /^</ && $line !~ /<br>$/) { # middle line
				$targets{$pathway}{Description} .= $line;
			} elsif ($line =~ /^<\/div>/) {
				$capture = undef;
	#            print "  Description: $targets{$pathway}{Description}\n";
			} else {
				print "Failed to comprehend line: $line\n";
			}
		} elsif ($capture eq 'Class') {
			($targets{$pathway}{Class}) = ($line =~ />([^<>]+)<br>$/);
			$capture = undef;
	#        print "  Class: $targets{$pathway}{Class}\n";
		} elsif ($capture eq 'Gene') {
			foreach $field (split /<\/tr><tr>/, $line) {
				next unless $field =~ /href/;   # gene entries have links
				my $id;
				if ($orgcode eq 'dme') {
					($id) = ($field =~ />(Dmel_CG\d+)<\/a>/);
				} else {
					($id) = ($field =~ />(\d+)<\/a>/);
				}
				my ($gene) = ($field =~ /<td align="left">([^\[\(]+)/);
				$gene =~ s/<a href="\/dbget-bin\/www_bget\?ko:K\d+">K\d+<\/a>//g;
				$gene =~ s/\s*$//;
				push @{ $targets{$pathway}{Genes} }, "$id ($gene)";
				push @{ $targets{$pathway}{FGenes} }, "$id\t$gene";
			}
			$capture = undef;
		} elsif ($line !~ /\S/) {
			$capture = undef;   # blank lines definitely end capture, especially for 'Genes'
		}
	}
}

open OUT, "> ${orgcode}_pathway_finaldata.txt";
open FLAT, "> ${orgcode}_pathway_flat_data.txt";
print FLAT "Pathway\tPathName\tEntrezGene\tGeneName\n";
foreach my $pathway (sort keys %targets) {
    my $string = join "\t", map { $targets{$pathway}{$_} } qw/ Name Class Description /;
    my $genes = join ' // ', @{ $targets{$pathway}{Genes} };
    print OUT "$pathway\t$string\t$genes\n";
    print FLAT "$pathway\t$targets{$pathway}{Name}\t$_\n" foreach @{ $targets{$pathway}{FGenes} };
}
close OUT;
close FLAT;


chomp(my $end = `date`);
print "NewKEGG.pl '$orgcode' complete: $end\n";
exit;



sub query_dump {
	my ($responseContent, $dumpFile) = @_;
	$page1 = $page2 = $responseContent;    # $page1, $page2 are global
	$page1 =~ s/<br>/<br>\n/mg;
	$page2 =~ s/\n//mg;
	open OUT, "> $dumpFile";
	print OUT $page1;
	close OUT;
}

