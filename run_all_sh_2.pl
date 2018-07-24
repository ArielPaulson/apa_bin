#p! /usr/bin/perl
####################
# run_all.pl - given targets.txt, run tophat or cufflinks
# By Madelaine Gogol / Ariel Paulson
# 7/2011
####################

#run_all.pl targets.txt

my $trim = 0;
my $multi = 1;
my $donovel = 0;
my $finaldir = "~/bioan/Yu/sachiko_yamanaka/MOLNG-40_20120110/data";

my ($targets, $mode) = @ARGV;   # mode = 1, 2, or 3: 1 = tophat; 2 = cufflinks; 3 = post-cufflinks; 4 = follow-ups (later merge with 1)
my $bowtie = "/n/site/inst/Linux-x86_64/bioinfo/bowtie/bowtie-0.12.7/bowtie";
my $tophat = "/n/site/inst/Linux-x86_64/bioinfo/tophat/tophat-1.4.0.Linux_x86_64/tophat";
my $cufflinks = "/n/site/inst/Linux-x86_64/bioinfo/cufflinks/cufflinks-1.3.0.Linux_x86_64/cufflinks";
my $fasta = "/n/projects/apa/stuff/bowtie_building/mm9/mm9.fa";
my $fastaidx = "/n/projects/apa/stuff/bowtie_building/mm9/mm9.fa.fai";
die "No mode specified!\n" unless $mode;

open ALL, "> launch.$mode.sh";
open(targets,"$targets");

while(<targets>)
{
    chomp;
    my ($path,$lane,$barcode,$alias,$genome,$bed,$gtf) = split("\t",$_);
    if($path ne "path")
    {

#	print "$. $alias\n";
	my $outputdir = "data/$alias";
	open OUT, "> $alias.$mode.sh";
	    
	if ($mode == 1) {
	    
	    print "mkdir -p $outputdir\n";
	    system "mkdir -p $outputdir";
	    if ($donovel) {
		print "mkdir -p $outputdir/novel\n";
		system "mkdir -p $outputdir/novel";
	    }
	    
	    my @PE = glob "$path/s_${lane}_*_$barcode.fastq.gz";

	    if ($PE[1]) {     # paired-end

		print "End 1: $PE[0]\nEnd 2: $PE[1]\n";   # hopefully!
		
		my $suffix = "txt.gz";
		my $suffix = "txt";

		if ($trim) {
#		    print OUT "gunzip -c $PE[0] | fastx_trimmer -Q 33 -l $trim | gzip -c > data/${alias}_1.txt.gz\n"; 
#		    print OUT "gunzip -c $PE[1] | fastx_trimmer -Q 33 -l $trim | gzip -c > data/${alias}_2.txt.gz\n";
		    print OUT "gunzip -c $PE[0] | fastx_trimmer -Q 33 -l $trim > data/${alias}_1.txt\n"; 
		    print OUT "gunzip -c $PE[1] | fastx_trimmer -Q 33 -l $trim > data/${alias}_2.txt\n";
		} else {
		    print OUT "gunzip -c $PE[0] > data/${alias}_1.txt\n"; 
		    print OUT "gunzip -c $PE[1] > data/${alias}_2.txt\n";
		}

		my $bowtieparams = "-p 2 -m $multi -k $multi-o $outputdir --mate-inner-dist 200 --mate-std-dev 70 --segment-length 35 --segment-mismatches 2 $genome";
		my $com = "$bowtie $bowtie_params data/${alias}_1.$suffix data/${alias}_2.$suffix && mail -s $alias.tophat apa\@stowers.org < /dev/null";
		print OUT "$com\n"; # system $com;
	    
		my $tophat_params = "-G $gtf -p 2 -g $multi -o $outputdir --mate-inner-dist 200 --mate-std-dev 70 --segment-length 35 --segment-mismatches 2 $genome";
		my $com = "$tophat $tophat_params data/${alias}_1.$suffix data/${alias}_2.$suffix";
		print OUT "$com\n"; # system $com;
		
	    } else {      # single-end
		
		print "End 1: $PE[0]\n";   # hopefully!
		
		my $suffix = "txt.gz";
		my $suffix = "txt";

		if ($trim) {
#		    print OUT "gunzip -c $PE[0] | fastx_trimmer -Q 33 -l $trim | gzip -c > data/${alias}_1.txt.gz\n"; 
		    print OUT "gunzip -c $PE[0] | fastx_trimmer -Q 33 -l $trim > data/${alias}_1.txt\n"; 
		} else {
		    print OUT "gunzip -c $PE[0] > data/${alias}_1.txt\n"; 
		}
		
		my $bowtieparams = "-p 2 -m $multi -k $multi -o $outputdir $genome";
		my $com = "$bowtie $bowtie_params data/${alias}_1.$suffix && mail -s $alias.tophat apa\@stowers.org < /dev/null";
#		print OUT "$com\n"; # system $com;
		
		my $tophat_params = "-G $gtf -p 2 -g $multi -o $outputdir --segment-length 35 --segment-mismatches 2 $genome";
		my $com = "$tophat $tophat_params data/${alias}_1.$suffix";
		print OUT "$com\n"; # system $com;
	    }
	    
	    #rename bam, index
	    my $com = "mv $outputdir/accepted_hits.bam $outputdir/$alias.bam\n";
	    print OUT "$com\n"; # system $com;
	    my $com = "samtools index $outputdir/$alias.bam\n";
	    print OUT "$com\n"; # system $com;
	    
	    my $com = "mail -s '$alias stage1 complete' apa\@stowers.org < /dev/null";
	    print OUT "$com\n"; # system $com;
	    
	} elsif ($mode == 2) {
	    
	    my $com = "perl uxonic_coverage.pl /n/data1/genomes/bowtie-index/mm9/Ens_63/mm9.Ens_63.uxons.bed $outputdir/$alias.bam";  # generates gene_uxonic_coverage.txt
	    print OUT "$com\n"; # system $com;
	    
	    my $cuffcall = $multi > 1 ? "$cufflinks -p 2 -u" : "$cufflinks -p 2";

#	    my $srcdir = '~/bioan/Yu/sachiko_yamanaka/MOLNG-40_20120110/data';
	    my $com = "$cuffcall -o $outputdir -G $gtf --max-mle-iterations 10000 $outputdir/$alias.bam";
	    print OUT "$com\n"; # system $com;
	    
	    if ($donovel) {
		my $com = "$cuffcall -o $outputdir --max-mle-iterations 10000 $outputdir/$alias.bam";
		print OUT "$com\n"; # system $com;
		
		my $com = "perl novel_cufflinks_genes.pl mm9.Ens_63.gff $outputdir/novel/transcripts.gtf";  # generates novel_transcripts.gtf
		print OUT "$com\n"; # system $com;
	    }
	    
	    my $com = "mail -s '$alias stage2 complete' apa\@stowers.org < /dev/null";
	    print OUT "$com\n"; # system $com;
	    
	} elsif ($mode == 3) {
	    
	    my $com = "java -jar /n/site/inst/Linux-x86_64/bioinfo/picard/current/CollectAlignmentSummaryMetrics.jar VALIDATION_STRINGENCY=LENIENT REFERENCE_SEQUENCE=$fasta INPUT=$outputdir/$alias.bam OUTPUT=$outputdir/picard.txt\n";
	    print OUT "$com\n"; # system $com;

	    #bw
	    my $com = "genomeCoverageBed -split -bg -ibam $outputdir/$alias.bam -g $fastaidx > $outputdir/$alias.bg\n"; 
	    print OUT "$com\n"; # system $com;
	    my $com = "bedGraphToBigWig $outputdir/$alias.bg $fastaidx $outputdir/$alias.bw\n";	
	    print OUT "$com\n"; # system $com;

	    #copy relevant files to appropriate locations
	    my $com = "cp $outputdir/picard.txt $finaldir/$alias.picard.txt\n";
	    print OUT "$com\n"; # system $com;
	    my $com = "cp $outputdir/gene_uxonic_coverage.txt $finaldir/$alias.gene_uxonic_coverage.txt\n";
	    print OUT "$com\n"; # system $com;
	    my $com = "cp $outputdir/$alias.bw $finaldir\n";
	    print OUT "$com\n"; # system $com;
	    my $com = "cp $outputdir/$alias.bam $finaldir\n";
	    print OUT "$com\n"; # system $com;
	    my $com = "cp $outputdir/$alias.bam.bai $finaldir\n";
	    print OUT "$com\n"; # system $com;
	    my $com = "cp $outputdir/unmapped_left.fq.z $finaldir/$alias.unmapped_left.fq.z\n";
	    print OUT "$com\n"; # system $com;
	    my $com = "cp $outputdir/genes.fpkm_tracking $finaldir/$alias.genes.fpkm_tracking\n";
	    print OUT "$com\n"; # system $com;
	    my $com = "cp $outputdir/isoforms.fpkm_tracking $finaldir/$alias.isoforms.fpkm_tracking\n";
	    print OUT "$com\n"; # system $com;
	    my $com = "cp $outputdir/transcripts.gtf $finaldir/$alias.transcripts.gtf\n";
	    print OUT "$com\n"; # system $com;

	    if ($donovel) {
		my $com = "mkdir -p $finaldir/novel\n";
		print OUT "$com\n"; # system $com;
		my $com = "cp $outputdir/novel/genes.fpkm_tracking $finaldir/novel/$alias.genes.fpkm_tracking\n";
		print OUT "$com\n"; # system $com;
		my $com = "cp $outputdir/novel/isoforms.fpkm_tracking $finaldir/novel/$alias.isoforms.fpkm_tracking\n";
		print OUT "$com\n"; # system $com;
		my $com = "cp $outputdir/novel/transcripts.gtf $finaldir/novel/$alias.transcripts.gtf\n";
		print OUT "$com\n"; # system $com;
		my $com = "cp $outputdir/novel/novel_transcripts.gff $finaldir/novel/$alias.novel_transcripts.gtf\n";
		print OUT "$com\n"; # system $com;
	    }

	    my $com = "mail -s '$alias stage3 complete' apa\@stowers.org < /dev/null";
	    print OUT "$com\n"; # system $com;
	    
	} elsif ($mode == 4) {
	    
	    if (0) {
	    ## get multireads + unmappables
	    my @PE = glob "$path/s_${lane}_*_$barcode.fastq.gz";
	    print "End 1: $PE[0]\nEnd 2: $PE[1]\n";   # hopefully!
	    
	    print "gunzip -c $PE[0] | fastx_trimmer -Q 33 -l 70 | gzip -c > data/s_${lane}_1.txt.gz\n"; 
	    my $e1 = `gunzip -c $PE[0] | fastx_trimmer -Q 33 -l 70 | gzip -c > data/s_${lane}_1.txt.gz`;
	    print "E1 done: $e1\n";
	    print "gunzip -c $PE[1] | fastx_trimmer -Q 33 -l 70 | gzip -c > data/s_${lane}_2.txt.gz\n"; 
	    my $e2 = `gunzip -c $PE[1] | fastx_trimmer -Q 33 -l 70 | gzip -c > data/s_${lane}_2.txt.gz`;
	    print "E2 done: $e2\n";
	    sleep 1;
	    
#	    my $params = "-G $gtf -p 4 -g 1 -o $outputdir --mate-inner-dist 200 --mate-std-dev 70 --segment-length 35 --segment-mismatches 2 $genome";
#	    my $com = "$tophat $params data/${alias}_1.txt.gz data/${alias}_2.txt.gz"; 
	    
	    my $params = "";
	    my $com = "bowtie ";

	    print OUT "$com\n"; # system $com;
	    
	    ## EXA
	    print "gunzip -c $path/s_${lane}_1_$barcode.export.txt.gz > data/$alias.export.1.txt\n"; 
	    my $e1 = `gunzip -c $path/s_${lane}_1_$barcode.export.txt.gz > data/$alias.export.1.txt`; 
	    print "E1 done: $e1\n";
	    print "gunzip -c $path/s_${lane}_2_$barcode.export.txt.gz > data/$alias.export.2.txt\n"; 
	    my $e2 = `gunzip -c $path/s_${lane}_2_$barcode.export.txt.gz > data/$alias.export.2.txt`; 
	    print "E2 done: $e2\n";
	    sleep 1;
	    
#	    system "perl Eland_Export_Analysis.pl data/
	    }
	    
	    
	    ## summarize alignments
	    my $com = "perl BamSummary.pl $outputdir/$alias.bam";
	    print OUT "$com\n"; # system $com;

	    my $com = "mail -s '$alias stage4 complete' apa\@stowers.org < /dev/null";
	    print OUT "$com\n"; # system $com;
	}
	close OUT;
	print ALL "./$alias.$mode.sh &\n";
    }
}
close ALL;
system "chmod a+x *.sh";
exit;
