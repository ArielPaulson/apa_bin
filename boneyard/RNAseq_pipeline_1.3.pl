#! /usr/bin/perl
use IO::CaptureOutput qw / capture_exec /;
use Getopt::Long;
use apa_routines;
use strict;
use Cwd;

####################
# RNAseq_pipeline.pl // given targets.txt, run bowtie, tophat, cufflinks, and everything else // can specify how many lanes to run concurrently
# analysis code by Madelaine Gogol // pipeline by APA
# 8/2010
####################

# RNAseq_pipeline.pl targets.txt num_concurrent_procs

##################### USER EDIT

my $run_bowtie = 0;						# 1/0: toggle bowtie calls (comparison purposes only)
my $run_htseq = 0;						# 1/0: toggle htseq calls (comparison purposes only)

my %static;
$static{spikescript} = '/n/projects/apa/stuff/mcm_spikes.R';		# R script for analyzing spikes

my $coverageBed ='/n/local/stage/bedtools/bedtools2-2.19.0/bin/coverageBed';

###################### END EDIT

my ($targets, $restart_at_cufflinks);
my ($concurrent, $nthreads) = (1, 1);	# defaults

GetOptions(
	"f=s" => \$targets, 
	"c=s" => \$concurrent, 
	"t=s" => \$nthreads, 
	"r=s" => \$restart_at_cufflinks
);

die "No targets file specified!\n" unless -e $targets;

my $base_dir = cwd();
my $logspacer = '********************';
my @fieldnames = qw/ file run_name output_dir multimatch spikes bowtie_index bed gff gtf chrom_sizes /;
my %filepos = map {($_=>1)} (0,5..9);		# positions (in @fieldnames) corresponding to files
my ($index_path, $count, @logs, @logsCL, %lost);
#my @bowtie_components = ('.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt');
my @bowtie_components = ('.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt', '.fa.fai');


##### TEST STATIC FILES

print "\nStatic files:\n";
foreach my $key (keys %static) {
	if (-e $static{$key}) {
		print "Validated: $static{$key}\n";
	} else {
		print "** Failed: $static{$key}\n";
#		$lost{$static{$key}} = 1;
	}
}

##### READ TARGETS FILE

open targets, "$targets";
while (<targets>) {
	
	$_ =~ s/[\n\r\"]//g;
	$count++;
	next if $count == 1;		# skip header
	my @data = split "\t", $_;	# $file $run_name $output_dir $spikes $bowtie_idx $bed $gff $gtf $chrom_sizes
	my $run_name = @data[1];
	print "\n$run_name:\n";
	my $run_dir;
	
	##### TEST FILES / ADJUST PATHS / CREATE RUN DIRECTORIES (per lane)
	
	foreach my $i (0..$#data) {
		if ($filepos{$i}) {		# file; test it
			my ($oldfile, $newfile) = ($data[$i], '');
			if ($i == 5) {		# bowtie index
				$oldfile = "$base_dir/$oldfile" unless $oldfile =~ /^\//;	# if not fully pathed, append to current path
				$index_path = $1 if $oldfile =~ /^(.*)\//;			# greedy matching returns full path to bowtie index
				foreach my $suffix (@bowtie_components) {			# required bowtie index components
					my $component = $oldfile.$suffix;
					if (-e $component) {
						print "Validated: $component\n";
					} else {
						print "** Failed: $component\n";
						$lost{$component} = 1;
					}
				}
				$data[$i] = $oldfile;
			} else {
				if ($oldfile !~ /^\// && -e "$index_path/$oldfile") {		# if not fully pathed, test bowtie index path first
					$newfile = "$index_path/$oldfile";
					print "Validated: $newfile\n";
				} elsif ($oldfile !~ /^\// && -e "$base_dir/$oldfile") {	# if not in index path, test working dir
					$newfile = "$base_dir/$oldfile";
					print "Validated: $newfile\n";
				} elsif (-e $oldfile) {						# otherwise, test original path
					$newfile = $oldfile;
					print "Validated: $newfile\n";
				} else {							# otherwise, fail
					print "** Failed: $oldfile\n";
					$lost{$oldfile} = 1;
				}
				$data[$i] = $newfile;						# after tests, replace with repathed version
			}
		} elsif ($i == 2) {		# output directory; test it
			$data[$i] = "$base_dir/$data[$i]" unless $data[$i] =~ /^\//;	# if output_dir not fully pathed, append to current path
			$run_dir = "$data[$i]/$run_name";
			unless (-d $run_dir) {
				die "Cannot validate run directory '$run_dir' to work in: $!.  Halting.\n" if $restart_at_cufflinks;
				system "mkdir -p $run_dir";
				if (-d $run_dir) {	# did mkdir succeed?
					print "Created output directory '$run_dir'\n";
				} else {
					die "Failed to create output directory '$run_dir': $!.  Halting.\n";
				}
			}
		}
	}
	
	my $logfile = "$run_dir/$run_name.log";
	push @logs, $logfile;
	my $logfileCL = "$run_dir/${run_name}_cufflinks.log";
	push @logsCL, $logfileCL;
	
	##### INITIALIZE LOG FILES (per lane)
	
	open OUT, "> $logfile";
	print OUT "base_dir\t$base_dir\n";
	print OUT "run_dir\t$run_dir\n";
	print OUT "final_dir\t$data[2]/final\n";
	print OUT "$fieldnames[$_]\t$data[$_]\n" foreach (0..$#data);
	close OUT;
	
	open OUT, "> $logfileCL";
	close OUT;
}
close targets;

## if any files not validated, die
die "\n\nSome essential files could not be validated!  Exiting.\n\n" if %lost;

## figure out how many run blocks are needed for $runs total jobs @ $concurrent degree of parallelization
my $runs = scalar @logs;
my $blocks = int ($runs / $concurrent);
$blocks++ if ($runs % $concurrent > 0);

##### EXECUTE LANE PROCESSES IN BATCHES

my $already = 0;
my $child = 0;
foreach my $i (1..$blocks) {	# run blocks
	my @allpids;
	foreach my $j (1..$concurrent) {	# concurrent processes (per block)
		$already++;
		last if ($already > $runs);	# stop when done with all processes specified in targets.txt
		{	# begin fork retry block
			my $pid = fork();
			if ($pid) {		# parent; collect child PIDs
				push @allpids, $pid;
			} elsif ($pid == 0) {	# child; run lane-processing subroutine and exit
				&process_lane($already-1);
				exit;
			} else {		# error
				warn "Couldn't fork: $!";	# low resources, usually
				sleep 900;	# wait 15 minutes
				redo;		# then retry fork
			}
		}	# end fork retry block
	}
	&babysit($i, @allpids);		# parent will wait to execute next run block until all pending PIDs close out
}
exit;


###################################################################  SUBROUTINES  ###########################################################################


sub babysit {
	my ($I, @running) = @_;
	my @pidnames = @running;	# @running gets overwritten; @pidnames does not
	my $waiting;
	do {
		sleep 60;	# check pid status every minute
		$waiting = undef;
		foreach my $k (0..$#running) {
			next if $running[$k] == -1;			# ignore already-completed processes
			$running[$k] = waitpid $running[$k], 0;		# gets itself, or converts to -1
			$waiting = 1;					# flag: PIDs still running, so keep waiting
		}
	} while $waiting;
	my $signals = join "\n", map { "$pidnames[$_] = $running[$_]" } (0..$#running);
	my $now = timestamp('FULL');
	print "\n$logspacer Pipeline host PID $$ says: BATCH $I COMPLETE: $now\n$signals\n";
}


sub execute {
	my ($COM, $CAP, $LOG, $PROC) = @_;
	my $now = timestamp('FULL');
	&logreport("\n$logspacer $PROC | $now\n$COM\n", $LOG);	# reports command execution to logfile and screen
	if ($CAP) {
		my ($success, $messages) = capture_exec($COM);	# capture any STDOUT/STDERR from command and report this also
		if ($messages) {
			if ($CAP == 1) {
				&logreport("MESSAGES:\n$messages\n", $LOG);
			} else {
				open OUT, ">> $LOG";	# write to log but not screen (i.e. cufflinks output)
				print OUT "MESSAGES:\n$messages\n";
				close OUT;
			}
		}
	} else {
		system $COM;	# execute command with no STDOUT/STDERR capturing or reporting
	}
}


sub process_lane {	# the actual analysis code
	
	### Initialize
	
	my $instance = shift;
	my $logfile = $logs[$instance];
	my $logfileCL = $logsCL[$instance];
	my %vars;
	
	### Collect run parameters from logfile header and put in %vars

	open IN, $logfile or die "Cannot open logfile '$logfile'!\n";
	while (<IN>) {
		chomp;
		my ($var, $val) = split /\t/, $_;
		$vars{$var} = $val;
	}
	close IN;
	
	my $procID = "$vars{run_name} process PID $$";
	my @capture = (1, $logfile, $procID);		# static variable sets for &execute (with STDOUT/STDERR capture)
	my @captureCL = (2, $logfileCL, $procID);	# static variable sets for &execute (with STDOUT/STDERR capture to cufflinks logfile but NOT to screen!)
	my @nocapture = (0, $logfile, $procID);		# static variable sets for &execute (with no capture)
	
	my $now = timestamp('FULL');
	&logreport("\n$logspacer $procID\n$now: Initialized.\n", $logfile);
	
	### Set up some variables
	
	my $basedir = $vars{base_dir};
	my $outputdir = $vars{run_dir};
	my $bowtie_outputdir = "$outputdir/bowtie";
	my $tophat_outputdir = "$outputdir/tophat";
	my $cufflinks_outputdir = "$tophat_outputdir/cufflinksG";
	my $finaldir = $vars{final_dir};
	
	#my $bowtie = "bowtie";
	my $bowtie = "/n/site/inst/Linux-x86_64/bioinfo/bowtie/bowtie-0.12.5/bowtie";
	my $cufflinks = "/n/site/inst/Linux-x86_64/bioinfo/cufflinks-1.0.3.Linux_x86_64/cufflinks";
	my $tophat = "/n/site/inst/Linux-x86_64/bioinfo/tophat-1.3.1.Linux_x86_64/tophat";
	
	### Create subdirectories
	
	foreach my $dir ($tophat_outputdir, $cufflinks_outputdir, $finaldir) {
		&execute("mkdir -p $dir", @nocapture) unless -d $dir;
	}
	if ($run_bowtie) {
		&execute("mkdir -p $bowtie_outputdir", @nocapture) unless -d $bowtie_outputdir;
	}
	sleep 5;	# let concurrent processes get over this step together -- makes STDOUT viewing a little better!!
	
	### Run tophat, cufflinks

	unless ($restart_at_cufflinks) {
	    if ($vars{file} =~ /\.gz/) {
		&execute("gunzip -c $vars{file} | $tophat --GFF $vars{gff} -p $nthreads -g $vars{multimatch} -o $tophat_outputdir $vars{bowtie_index} -", @capture);
	    } else {
		&execute("$tophat --GFF $vars{gff} -p $nthreads -g $vars{multimatch} -o $tophat_outputdir $vars{bowtie_index} $vars{file}", @capture);
	    }
	}
	
	chdir $tophat_outputdir or die "Cannot chdir into tophat directory: $!\n";
	my $here1 = cwd();
	logreport ("Working in '$here1'\n", $logfile);
	&execute("sort -k 3,3 -k 4,4n accepted_hits.sam > sorted_accepted_hits.sam", @capture) unless $restart_at_cufflinks;

	chdir $cufflinks_outputdir or die "Cannot chdir into cufflinks directory: $!\n";
	my $here1 = cwd();
	logreport ("Working in '$here1'\n", $logfile);
	&execute("$cufflinks -u -g $vars{gtf} $tophat_outputdir/sorted_accepted_hits.sam", @captureCL);
	
	chdir $basedir or die "Cannot chdir into base directory: $!\n";
	my $here2 = cwd();
	logreport ("\nWorking in '$here2'\n", $logfile);
	
	### Run post-cufflinks stuff
	
	if ($vars{spikes}) {
		#count spikes
		&execute("cut -f 3 $tophat_outputdir/sorted_accepted_hits.sam | sort | uniq -c | sort -k 2 | perl -p -e \"s/^ +//g; s/ /\t/g;\" > $tophat_outputdir/chrom_count.txt", @capture); 
#		&execute("R --vanilla --slave -f $static{spikescript} --args $tophat_outputdir", @capture);
	}
	
	#tophat bam file
	&execute("samtools import $vars{bowtie_index}.fa.fai $tophat_outputdir/sorted_accepted_hits.sam $tophat_outputdir/accepted_hits.bam", @capture);
	&execute("samtools sort $tophat_outputdir/accepted_hits.bam $tophat_outputdir/accepted_hits.sorted", @capture);
	&execute("samtools index $tophat_outputdir/accepted_hits.sorted.bam", @capture);
	&execute("samtools flagstat $tophat_outputdir/accepted_hits.sorted.bam > $tophat_outputdir/accepted_hits.flagstat.txt", @capture);
	
	#tophat wig file
	&execute("perl -i -p -e 's/^(.*)/chr$1/ unless \$_ =~ /^(chr|track)/' $tophat_outputdir/coverage.wig", @capture);
	if ($vars{spikes}) {
		&execute("sed \"/^chr[DAP|LYS|PHE|THR|TRP].*/d\" $tophat_outputdir/coverage.wig | tail -n +2 > $tophat_outputdir/coverage_ucsc_nh.wig", @capture);
	} else {
		&execute("tail -n +2 $tophat_outputdir/coverage.wig > $tophat_outputdir/coverage_ucsc_nh.wig", @capture);
	}
	
	#bw
	&execute("wigToBigWig -clip $tophat_outputdir/coverage_ucsc_nh.wig $vars{chrom_sizes} $tophat_outputdir/$vars{run_name}.bw", @capture);	
	
	#tophat counts
	&execute("$coverageBed -abam $tophat_outputdir/accepted_hits.sorted.bam -b $vars{bed} | sort -r -n -k7 > $tophat_outputdir/counts.txt", @capture);
	&execute("sed -i -e '1ichrom	start	end	id	score	strand	count	numbasescov	genelength	fracgenecov' $tophat_outputdir/counts.txt", @capture);
	
	if ($run_htseq) {
		#htseq counts tophat
		&execute("python2.6 -m HTSeq.scripts.count -m intersection-nonempty --type=exon -s no $tophat_outputdir/sorted_accepted_hits.sam $vars{gtf} --quiet > $tophat_outputdir/htseq_counts.txt", @capture);
	}
		
	if ($run_bowtie) {
		#bowtie
		&execute("$bowtie $vars{bowtie_index} --max $bowtie_outputdir/multi.fq --un $bowtie_outputdir/unmapped.fq --sam -v 2 -k $vars{multimatch} -p $nthreads -m $vars{multimatch} $vars{file} > $bowtie_outputdir/sorted_accepted_hits.sam", @capture); 
		
		#bowtie bam file
		&execute("samtools import $vars{bowtie_index}.fa.fai $bowtie_outputdir/sorted_accepted_hits.sam $bowtie_outputdir/accepted_hits.bam", @capture);
		&execute("samtools sort $bowtie_outputdir/accepted_hits.bam $bowtie_outputdir/accepted_hits.sorted", @capture);
		&execute("samtools index $bowtie_outputdir/accepted_hits.sorted.bam", @capture);
		&execute("samtools flagstat $bowtie_outputdir/accepted_hits.sorted.bam > $bowtie_outputdir/accepted_hits.flagstat.txt", @capture);
		
		#bowtie counts
		&execute("$coverageBed -abam $bowtie_outputdir/accepted_hits.sorted.bam -b $vars{bed} | sort -r -n -k7 > $bowtie_outputdir/counts.txt", @capture);
		&execute("sed -i -e '1ichrom	start	end	id	score	strand	count	numbasescov	genelength	fracgenecov' $bowtie_outputdir/counts.txt", @capture);
		
		if ($run_htseq) {
			#htseq counts bowtie
			&execute("python2.6 -m HTSeq.scripts.count -m intersection-nonempty --type=exon -s no $bowtie_outputdir/sorted_accepted_hits.sam $vars{gtf} --quiet > $bowtie_outputdir/htseq_counts.txt", @capture);
		}
	}
	
	#copy relevant files to appropriate locations
	if ($vars{spikes}) {
#		&execute("cp $tophat_outputdir/spikes.pdf $finaldir/$vars{run_name}_spikes.pdf", @capture);
#		&execute("cp $tophat_outputdir/spike_counts.txt $finaldir/$vars{run_name}_spike_counts.txt", @capture);
	}
	if ($run_bowtie) {
		&execute("cp $bowtie_outputdir/counts.txt $finaldir/$vars{run_name}_bowtie_counts.txt", @capture);
	}
	&execute("cp $tophat_outputdir/$vars{run_name}.bw $finaldir", @capture);
	&execute("cp $tophat_outputdir/counts.txt $finaldir/$vars{run_name}_tophat_counts.txt", @capture);
	&execute("cp $tophat_outputdir/cufflinksG/genes.fpkm_tracking $finaldir/$vars{run_name}.genes.fpkm_tracking", @capture);
	&execute("cp $tophat_outputdir/cufflinksG/isoforms.fpkm_tracking $finaldir/$vars{run_name}.isoforms.fpkm_tracking", @capture);
	
	my $now = timestamp('FULL');
	&logreport("\n$now: Processing Complete!\n", $logfile);
}
