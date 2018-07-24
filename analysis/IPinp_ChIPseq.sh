#!/bin/bash

## Intended to be run once for a SINGLE input, and ONE OR MORE IPs related to that input.
## Must be run from a 'code' directory, i.e. /path/to/cbio.xxx.nnn/code

## IPs and inputs must be directory names in ../data/output/; each directory must contain a 'bowtie2.bam'
## All 'bowtie2.bam' files should have 3 accompanying files: bowtie2.bam.bai, bowtie2.idxstats.txt, bowtie2.RChipQC.RData
## If using pooled inputs or IPs, have separate pool dirs like 'Input_pooled' with pooled bam named 'Input_pooled/bowtie2.bam'; also have 3 accompanying files

inp=$1   # single input
IPs=$2   # single IP, or comma-delim strimg of IPs
geno=$3  # MACS2 'genome' param
cname=$4 # name (pathless prefix) for RChance output (will write to ../data/chance)

IPs=($(echo $IPs | sed 's/,/\n/g'))  # break up comma-delim string of IPs, if multiple IPs

code=$(pwd)
cbio=${code%/code}
out=$cbio/data/output
macs=$cbio/data/macs2
chance=$cbio/data/chance

if [ ! -d $macs ]; then mkdir $macs; fi
if [ ! -d $chance ]; then mkdir $chance; fi

IPlist=''
for IP in ${IPs[@]}
do
    IPlist="$IPlist $IP,$out/$IP/bowtie2.RChipQC.RData"
    macs2 callpeak -t $out/$IP/bowtie2.bam -c $out/$inp/bowtie2.bam -g $geno -n $IP.macs2 --outdir $macs/$IP
done

/home/apa/local/bin/RChance.R ../data/chance/$cname $IPlist input=$inp,../data/output/$inp/bowtie2.RChipQC.RData &

