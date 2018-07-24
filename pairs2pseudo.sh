#!/bin/bash

## FUTURE: Perl script for simpler argument handling; input N bams, get N pseudoreplicates

## Usage: pairs2pseudo outprefix bam1 bam2 [--clobber]

## Given a output prefix (preferably '/path/to/pseudo') and 2 sorted replicate bams, creates the following:
## $prefix-1.bam
## $prefix-2.bam
## and associated bam indexes, idxstats.txt files

prefix=$1
bam1=$2
bam2=$3
clobber=$4

## Test command line arguments
if [ -z $prefix ] || [ -z $bam1 ] || [ -z $bam2 ]; then
    echo "Bad command line: use 'pairs2pseudo <prefix> <bam1> <bam2> <clobber>'";
    exit 1
fi
if [ ! -e $bam1 ]; then
    echo "Bam file 1 '$bam1' does not exist!";
    exit 1
else
    bam1=$(readlink -f $bam1)
fi

if [ ! -e $bam2 ]; then
    echo "Bam file 2 '$bam2' does not exist!";
    exit 1
else
    bam2=$(readlink -f $bam2)
fi

if [ -z $clobber ]; then
    clobber=0
else
    if [ "$clobber" == "1" ]; then clobber=1; else clobber=0; fi
fi

## Variables
bamname1=${bam1%.bam}
bamname2=${bam2%.bam}
prefix=$(readlink -f $prefix)  # in case it began with "../" or "~/" or suchlike
pseu1=$prefix-1
pseu2=$prefix-2

if [ -e $pseu1.idxstats.txt ] && [ -e $pseu2.idxstats.txt ] && [ $clobber -ne 1 ]; then
	echo "Pseudoreplicates already created!  Skipping."
	exit 0
fi

## Make temp directory
tmp="pairs2pseudo.$$.tmp"
if [ -d $tmp ]; then
    echo "Temp directory '$tmp' already exists!";
    exit 1
fi
mkdir $tmp
if [ ! -d $tmp ]; then
    echo "Could not make temp directory '$tmp'!";
    exit 1
fi
echo "Temp directory: $tmp"

## Ensure bam indexes and idxstats.txt files exist
if [ ! -e $bam1.bai ]; then
	echo "Indexing $bam1..."
	samtools index $bam1
fi
if [ ! -e $bamname1.idxstats.txt ]; then
	samtools idxstats $bam1 > $bamname1.idxstats.txt
fi
if [ ! -e $bam2.bai ]; then
	echo "Indexing $bam2..."
	samtools index $bam2
fi
if [ ! -e $bamname2.idxstats.txt ]; then
	samtools idxstats $bam2 > $bamname2.idxstats.txt
fi

## How many reads per pseudoreplicate?
N1=$(paste -s -d+ <(cut -f3 $bamname1.idxstats.txt) | bc)
N2=$(paste -s -d+ <(cut -f3 $bamname2.idxstats.txt) | bc)
NH=$(echo "($N1+$N2)/2" | bc)
echo -e "Bam 1:  $N1\nBam 2:  $N2\nPseudo: $NH each"

## Create pseudoreps
echo "Creating pooled sam..."
samtools view $bam1 > $tmp/sam
samtools view $bam2 >> $tmp/sam
samtools view -H $bam1 > $tmp/ps1.sam
samtools view -H $bam2 > $tmp/ps2.sam

echo "Creating pseudoreplicates..."
cd $tmp
cat sam | shuf | split --lines=$NH -
cat xaa >> ps1.sam
cat xab >> ps2.sam
samtools view -bS ps1.sam > ps1.bam
samtools view -bS ps2.sam > ps2.bam
samtools sort -m 10G ps1.bam ps1.sorted
samtools sort -m 10G ps2.bam ps2.sorted
cd ..

mv -f $tmp/ps1.sorted.bam $pseu1.bam
mv -f $tmp/ps2.sorted.bam $pseu2.bam

## Transfer pseudoreps to final location; index and idxstat
echo "Indexing pseudoreplicates..."
samtools index $pseu1.bam
samtools index $pseu2.bam
samtools idxstats $pseu1.bam > $pseu1.idxstats.txt
samtools idxstats $pseu2.bam > $pseu2.idxstats.txt

## Exit
rm -rf $tmp
echo "pairs2pseudo $outpath complete!"
exit 0
