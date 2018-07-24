#!/bin/bash

## Must be run from a 'code' directory, i.e. /path/to/cbio.xxx.nnn/code
## Not intended to be run at-once.  Provides generic ChIPseq workflow and commands.

for subdir in fastq output bigWig; do mkdir ../data/$subdir; done

cp ~/lbin/analysis/run_RNAseq.sh run.sh

## Core params
lab=""
usr=""
order=""
flowcell=""
geno=""
anno=""
direc=1
cores=1
mem="1G"
nomerge=0  # 0 to merge multiple lane fastqs per sample; 1 to process independently ('0' produces 'merge.sh' which must be run first) 

## RNAseq-specific params

~/lbin/sampleReportToLaunch $lab $usr $order $flowcell $geno $anno $direc $cores $mem $nomerge  # creates 'launch.sh' and perhaps 'merge.sh'


./merge.sh

./launch.sh




~/lbin/mergeAlignSummary --tophat ../data/output/*/align_summary.txt > all.align_summary.txt 
perl -i -pe 's!../data/output/!!' all.align_summary.txt 
perl -i -pe 's!/align_summary.txt!!' all.align_summary.txt 

#~/lbin/mergeRnaSeqMetrics -o all.RnaSeqMetrics.txt --dirname ../data/output/*/RnaSeqMetrics.txt 

~/lbin/mergeRnaSeqMetrics -o all.RnaSeqMetricsPlus.txt --dirname ../data/output/*/RnaSeqMetricsPlus.txt 

~/lbin/mergeMatrix -o all.htseq_counts.txt -k 1 -v 2 -h 0 -kn GeneID --dirname ../data/output/*/htseq.counts.txt

~/lbin/mergeMatrix -o all.cufflinks_fpkms.txt -k 1 -v 10 -h 1 -kn GeneID --dirname ../data/output/*/genes.fpkm_tracking

~/lbin/bamPosCount ../data/output/*/accepted_hits_primary.bam > all.unq_pos.txt


## Unaligned reads analysis

awk '{ if ($12>0.4) print $0 }' all.align_summary.txt > unaligned.txt

while read line
do
    samp=$(echo $line | cut -f1 -d' ')
    out=../data/output/$samp
    echo $out
    samtools view $out/accepted_hits_primary.bam | cut -f1 > $out/aln.idx
    ~/lbin/readFilter -i ../data/fastq/$samp.fastq.gz -o $out/unaligned.fastq.gz -s $out/aln.idx --exclude --destructive --first-field
    ~/lbin/readSample -f $out/unaligned.fastq.gz -o $out/unaligned.top10k.fastq.gz -n 10000
    ~/lbin/readSampleEval -f $out/unaligned.top10k.fastq.gz -oc $out/unaligned.top10k.eval -r Mus_musculus --keep
done < <(tail -n +2 unaligned.txt)
