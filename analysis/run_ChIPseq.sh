#!/bin/bash

## Must be run from a 'code' directory, i.e. /path/to/cbio.xxx.nnn/code
## NO PROVISIONS FOR PAIRED-END AT THE MOMENT -- IGNORES $FASTQ2 ARGUMENT

sample=$1  # sample name; becomes output dir under ../data/output/
fastq1=$2  # end 1 fastq, or comma-separated fastq list
fastq2=$3  # end 2 fastq(s) OR 'NA' IF SINGLE-END *******
geno=$4    # e.g. 'dm3'
anno=$5    # e.g. 'Ens_78'
direc=$6   # 1|0 if directional libraries -- ignored for ChIPseq ops -- see below
cores=$7   # N cores for bowtie2, samtools sort; default 1
mem=$8     # mem for samtools sort, e.g. "10G"; default 1G
run_star=$9

## As of the moment (2.2.0), bowtie2 cannot treat directional libraries differently.
## The '$direc' argument exists for compatibility with ~/lbin/sampleReportToLaunch, 
##  but is otherwise incompatible with ChIPseq operations -- should always be '0'
direc=0

if [ -z $cores ]; then cores=1; fi
if [ -z $mem ]; then mem='1G'; fi

code=$(pwd)
cbio=${code%/code}
outdir=$cbio/data/output/$sample
bigwig=$cbio/bigWig

if [ ! -d $outdir ]; then mkdir -p $outdir; fi
if [ ! -d $bigwig ]; then mkdir -p $bigwig; fi

## genome data path
gdir=/n/data1/genomes/bowtie-index/$geno
## annotation prefix
annots=$gdir/$anno/$geno.$anno

if [ "$run_star" == "0" ]; then
    
    ## See Tophat manual at https://ccb.jhu.edu/software/tophat/manual.shtml
    
    ## Bowtie2 genome index
    bti=$gdir/$geno
    ## Tophat bam prefix
    all=$outdir/bowtie2
    
else
    
    ## See STAR manual at https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
    
    ## STAR index
    stari=$gdir/$anno/STAR_${run_star}bp
    ## STAR bam prefix
    all=$outdir/Aligned.out
    
fi



if [ "$run_star" == "0" ]; then
    
    ##### Align reads (Bowtie2)

    ## 0: Make output directory
    mkdir -p $outdir
    ## 1: Bowtie2 alignment
    if [ "$fastq2" != "NA" ]; then fastq="-1 $fastq1 -2 $fastq2"; else fastq="-U $fastq1"; fi
    bowtie2 -p $cores --un-gz $all.unaligned.fastq.gz -x $bti $fastq -S $all.sam 2> $all.align_summary.txt
    
else
    
    ##### Align reads (STAR)
   
    ## 1: STAR alignment
    cd $outdir
    if [ "$fastq2" != "NA" ]; then fastq="$fastq1 $fastq2"; else fastq=$fastq1; fi
    STAR --genomeDir $stari --readFilesIn $fastq --readFilesCommand zcat --runThreadN $cores
    cd $code
    
fi
    
## 2: SAM -> BAM
samtools view -h -F 4 -bS $all.sam | samtools sort -@ $cores -m $mem -o $all.bam -
## 3: Index BAM file
samtools index $all.bam
## 4: Index stats
samtools idxstats $all.bam > $all.idxstats.txt
## 5. BedGraph -> bigWig conversion
/home/apa/local/bin/bam2track -b $all.bam -g $geno -n APM --BW
ln -sf $(readlink -f $all.bam) $bigwig/$sample.APM.bw
## 6. RChipQC.R
/home/apa/local/bin/RChipQC.R $all.bam 3 TRUE
## 7. Fastq reads, if single fastq given
if [ ! $fastq =~ "," ]; then /home/apa/local/bin/readCount $fastq --unique > $outdir/fqReads.txt; fi

rm -f $all.sam

## Email notification when complete
usr=$(whoami)
#mail -s "$sample complete" $usr@stowers.org < /dev/null

