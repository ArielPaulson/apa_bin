#!/bin/bash

sample=$1

geno=mm10
anno=Ens_87
cores=60
ram=100G

out=../data/output/$sample
mkdir $out


## Trim Reads
if [[ $sample =~ "Bioo_10ng" ]]; then
	~/lbin/trimReads -fq1 ../data/fastq/$sample.fastq.gz -c $cores -m $ram -t ~/cbio.cry.103/code/NextFlex.fa -x 4,4 -o $out/trimReads -min 15 -max 45
else
	~/lbin/trimReads -fq1 ../data/fastq/$sample.fastq.gz -c $cores -m $ram -o $out/trimReads -min 15 -max 45
fi


## Ribosomes first
~/lbin/alignToRibosomes -g $geno -a $anno -o $out/alignToRibosomes -c $cores -m $ram -fq1 $out/trimReads/trimmed.passing.fastq.gz --rl-split
## bp-level quantitation calls for later
pref=$out/alignToRibosomes/alignToRibosomes; for bam in $pref.RGL/*bp.bam; do echo "~/lbin/bamQuantitate -i $bam -o $bam.quant.txt --transcriptome" >> $pref.RGLQ.sh; done


## Functional RNA second
~/lbin/alignToFuncRNA -g $geno -a $anno -o $out/alignToFuncRNA -c $cores -m $ram -fq $out/alignToRibosomes/alignToRibosomes.unaligned.fastq.gz --rl-split
## bp-level quantitation calls for later
pref=$out/alignToFuncRNA/alignToFuncRNA; for bam in $pref.RGL/*bp.bam; do echo "~/lbin/bamQuantitate -i $bam -o $bam.quant.txt --transcriptome" >> $pref.RGLQ.sh; done


## Whole transcriptome third
~/lbin/alignToGenome -g $geno -a $anno -o $out/alignToGenome1T -c $cores -m $ram -aln ShortStack -fq1 $out/alignToFuncRNA/alignToFuncRNA.unaligned.fastq.gz --rl-split --smallRNA --transcriptome
## bp-level quantitation calls for later
pref=$out/alignToGenome1T/alignToGenome; for bam in $pref.RGL/*bp.bam; do echo "~/lbin/bamQuantitate -i $bam -o $bam.quant.txt -g mm10 -a Ens_87 --transcriptome" >> $pref.RGLQ.sh; done


## Rest of genome fourth
~/lbin/alignToGenome -g $geno -a $anno -o $out/alignToGenome2 -c $cores -m $ram -aln ShortStack -fq1 $out/alignToGenome1T/alignToGenome.unaligned.fastq.gz --rl-split --smallRNA
## bp-level quantitation calls for later
pref=$out/alignToGenome2/alignToGenome; for bam in $pref.RGL/*bp.bam; do echo "~/lbin/bamQuantitate -i $bam -o $bam.quant.txt -g mm10 -a Ens_87 --3pass" >> $pref.RGLQ.sh; done

