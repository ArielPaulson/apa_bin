#!/usr/bin/perl
use Cwd;

my $cores = 24;
my @K = (25,50,75,100,125);
my $base = cwd();
my $correct = $ARGV[0] ? 1 : 0;

open IN, 'targets.txt' or die;
while (<IN>) {
    $_ =~ s/[\n\r]+$//;
    my ($path, $lane, $alias) = split /\t/, $_;
    next if $path eq 'path';  # header
    next unless $alias =~ /^inornata/;  # testing only inornata libraries right now
    push @aliases, $alias;
    $data{$alias} = [$path, $lane];
    $files{$alias} = ["s_${lane}_1_NoIndex.fastq.gz", "s_${lane}_2_NoIndex.fastq.gz"];
#    $files{$alias} = [glob "$path/s_${lane}_*fastq.gz"];
#    print "$alias: @{ $files{$alias} }\n";
}
close IN;

foreach my $alias (@aliases) {
    foreach my $k (@K) {
	my $tag = "$alias.k$k";
	
	my $sdn;
	if ($k <= 31) {
	    $sdn = "$base/SOAPdenovo/SOAPdenovo31mer";
	} elsif ($k <= 63) {
	    $sdn = "$base/SOAPdenovo/SOAPdenovo63mer";
	} elsif ($k <= 127) {
	    $sdn = "$base/SOAPdenovo/SOAPdenovo127mer";
	}
#	print "$alias: $k: $sdn\n";

	mkdir "$base/$alias/k$k";
	chdir "$base/$alias/k$k";
	system "ln -s $data{$alias}->[0] $alias" unless -e $alias; # symlinked to get $alias to appear as library name
	open CFG, "> $tag.config";
	print CFG "max_rd_len=100\n";
	print CFG "[LIB]\n";
	print CFG "avg_ins=300\n";
	print CFG "reverse_seq=0\n";
	print CFG "asm_flags=3\n";
	print CFG "rank=1\n";
	print CFG "q1=./$alias/$files{$alias}->[0]\n"; # symlinked to get $alias to appear as library name
	print CFG "q2=./$alias/$files{$alias}->[1]\n";
	close CFG;
	
	system "date > $tag.runtimes";
	system "echo pregraph >> $tag.runtimes";
	system "$sdn pregraph -K $k -p $cores -R -s $tag.config -o $tag";
	system "date >> $tag.runtimes";
	system "echo contig >> $tag.runtimes";
	system "$sdn contig -g $tag";
	system "date >> $tag.runtimes";
	system "echo map >> $tag.runtimes";
	system "$sdn map -p $cores -s $tag.config -g $tag";
	system "date >> $tag.runtimes";
	system "echo scaff >> $tag.runtimes";
	system "$sdn scaff -p $cores -g $tag";
	system "date >> $tag.runtimes";
	
#	system "date >> $tag.runtimes";
#	system "echo GapCloser >> $tag.runtimes";
#	system "$base/Utils/GapCloser -b $tag.config -a $tag.scaff -t $cores -o $tag.GapClosed";
#	system "date >> $tag.runtimes";
	
	system "grep '>' $tag.contig > $tag.contig.headers";
	system "cut -d' ' -f3 $tag.contig.headers > $tag.contig.lengths";
	
    }
}

