#!/bin/bash
set -e

## Must be run from a 'code' directory, i.e. /path/to/cbio.xxx.nnn/code

sample=$1  # sample name; becomes output dir under ../data/output/
fastq1=$2  # end 1 fastq, or comma-separated fastq list
fastq2=$3  # end 2 fastq(s) OR 'NA' IF SINGLE-END *******
geno=$4    # e.g. 'dm3'
anno=$5    # e.g. 'Ens_78'
direc=$6   # 1|0 if directional libraries **** ASSUMES TRUSEQ DIRECTIONAL -- see Tophat flags a few lines below
cores=$7   # N cores for bowtie2, samtools sort; default 1
mem=$8     # mem for samtools sort, e.g. "10G"; default 1G
run_star=$9  # if running STAR not Tophat, give index size (i.e. 51, 76, 101, 151; that index must already be built)

if [ -z $direc ]; then direc=0; else direc=1; fi
if [ -z $cores ]; then cores=1; fi
if [ -z $mem ]; then mem='1G'; fi
if [ -z $run_star ]; then run_star=0; fi

if [ "$direc" == "1" ]; then
    tophat_flags="--library-type fr-firststrand"
    star_flags=""
else
    tophat_flags=""
    star_flags=""
fi

code=$(pwd)
cbio=${code%/code}
outdir=$cbio/data/output/$sample
bigwig=$cbio/data/bigWig

if [ ! -d "$outdir" ]; then mkdir -p $outdir; fi
if [ ! -d "$bigwig" ]; then mkdir -p $bigwig; fi
if [ "$fastq2" != "NA" ]; then fastq="$fastq1 $fastq2"; else fastq=$fastq1; fi

## genome data path
gdir=/n/data1/genomes/bowtie-index/$geno
## annotation prefix
annots=$gdir/$anno/$geno.$anno

if [ "$run_star" == "0" ]; then
    
    ## See Tophat manual at https://ccb.jhu.edu/software/tophat/manual.shtml
    
    ## Bowtie2 genome index
    bti=$gdir/$geno
    ## Bowtie2 GTF index
    gtfi=$annots.cuff.gtf.index/$geno.$anno.cuff
    ## Tophat bam prefix
    all=$outdir/accepted_hits
    ## primary bam
    pri=$outdir/accepted_hits_primary
    
else
    
    ## See STAR manual at https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
    
    ## STAR index
    stari=$gdir/$anno/STAR_${run_star}bp
    ## STAR bam prefix
    all=$outdir/Aligned.out
    ## primary bam
    pri=$outdir/Aligned.out.primary
    
fi







if [ "$run_star" == "0" ]; then
    
    ##### Align reads (Tophat)
    
    ## 1: Tophat alignment
    tophat -p $cores --transcriptome-index $gfti --no-novel-juncs --no-coverage-search -o $outdir $tophat_flags $bti $fastq
    
else
    
    ##### Align reads (STAR)
    
    ## 1: STAR alignment
    cd $outdir
    STAR --genomeDir $stari --readFilesIn $fastq --readFilesCommand zcat --runThreadN $cores --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --quantMode GeneCounts
    cd $code
    ## 1.5: STAR BAM sorting
    samtools view -h -F 4 -bS $all.bam | samtools sort -@ $cores -m $mem -o $outdir/sorted.bam -
    mv -f $outdir/sorted.bam $all.bam
    
fi



##### Post-alignment

## 2: Index BAM file
samtools index $all.bam
## 3: Extract primary alignments (only 1 align per read)
samtools view -h -F 256 $all.bam | samtools view -bS - > $pri.bam
## 4: Index primary-only BAM
samtools index $pri.bam


##### Split BAMs, if directional



##### QC

## 5: Count unique aligned positions
samtools view $pri.bam | cut -f3,4 | uniq > $pri.unq_pos.bed
## 6a. Count reads aligned to each chromosome (the fast way)
samtools idxstats $all.bam > $all.idxstats.txt
samtools idxstats $pri.bam > $pri.idxstats.txt


if [ "$run_star" == "0" ]; then
    
    ## 7. Picard CollectRnaSeqMetrics (primary alignments)
    ~/lbin/CollectRnaSeqMetricsPlus -g $geno -a $anno -b $pri.bam -m $mem -o $outdir/RnaSeqMetrics.txt
    
    ##### Quantitate gene models
    
    ## 8. Cufflinks on full bam (not primary alignments)
    cufflinks -p 4 -u -G $annots.cuff.gtf -o $outdir $all.bam --max-bundle-frags 1000000
    ## 9. HT-Seq count
    samtools view $pri.bam | htseq-count -s no -a 0 -m intersection-nonempty - $annots.cuff.gtf > $outdir/htseq.counts.txt
    # IF COUNTING MULTIREADS:
    # samtools view $pri.bam | perl -pe 's/\s+NH:i:\d+//' | htseq-count -s no -a 0 -m intersection-nonempty - $annots.cuff.gtf > $outdir/htseq.counts.txt
    
else
    
    ## 7. Picard CollectRnaSeqMetrics (primary alignments)
    ~/lbin/CollectRnaSeqMetricsPlus -g $geno -a $anno -b $pri.bam -m $mem -o $outdir/RnaSeqMetrics.txt -s samtools-sorted
    
fi


##### Visualize

if [ $direc -eq 1 ]; then strflag='--stranded --switch'; else strflag=''; fi
/home/apa/local/bin/bam2track -b $pri.bam -g $geno -n APM --BW $strflag
ln -sf $(readlink -f $pri.APM.bw) $bigwig/$sample.RPM.bw

~/lbin/readCount --uhisto $fastq1 > $fastq1.fqr &

## Email notification when complete
usr=$(whoami)
mail -s "$outdir complete" $usr@stowers.org < /dev/null

