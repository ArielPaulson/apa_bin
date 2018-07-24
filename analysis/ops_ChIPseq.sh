#!/bin/bash

## Must be run from a 'code' directory, i.e. /path/to/cbio.xxx.nnn/code
## Not intended to be run at-once.  Provides generic ChIPseq workflow and commands.

for subdir in fastq output bigWig chance; do mkdir ../data/$subdir; done

cp ~/lbin/analysis/run_ChIPseq.sh run.sh

## Core params
lab=""
usr=""
order=""
flowcell=""
geno=""
anno=""
direc=0
cores=1
mem="1G"
nomerge=0  # 0 to merge multiple lane fastqs per sample; 1 to process independently ('0' produces 'merge.sh' which must be run first) 

## ChIPseq-specific params
macs2geno=""
flankbp=500

 | sort -k 8nr,8nr | head -n 100000
~/lbin/sampleReportToLaunch $lab $usr $order $flowcell $geno $anno $direc $cores $mem $nomerge  # creates 'launch.sh' and perhaps 'merge.sh'

./merge.sh

./launch.sh

### NEED TO RUN MACS2 / create IP-input file etc..

mergedir=merged_peaks_$flankbp
~/lbin/mergePeaks -o $mergedir -f $flankbp -g $genome --max --strip .macs2_peaks ../data/macs2/*/*.narrowPeak

######## GENERALIZE THE BELOW -- RUN OFF IP COLUMN OF IP-INPUT FILE
ls -1 ../data/output | awk '{ print "IP="$1",../data/output/"$1"/bowtie2.RChipQC.RData" }'
~/lbin/ChIPBinCorr.R genome.1kbins 2677CtrlS2,../data/output/2677CtrlS2/bowtie2.RChipQC.RData 2677MTPD,../data/output/2677MTPD/bowtie2.RChipQC.RData 2677WTFull,../data/output/2677WTFull/bowtie2.RChipQC.RData HAMTPD,../data/output/HAMTPD/bowtie2.RChipQC.RData HAWTFull,../data/output/HAWTFull/bowtie2.RChipQC.RData InputCtrlS2,../data/output/InputCtrlS2/bowtie2.RChipQC.RData InputMTPD,../data/output/InputMTPD/bowtie2.RChipQC.RData InputWTFull,../data/output/InputWTFull/bowtie2.RChipQC.RData 

~/lbin/getRandomCoords -o random.bed -b $mergedir/final.bed -c ~/bwti/$geno/$geno.chrom.sizes -g ~/bwti/$geno/$geno.N-blocks.bed -s 10 --avoid
## OR ##
~/lbin/getRandomCoords -o random.bed -N 50000 -L 1000 -c ~/bwti/$geno/$geno.chrom.sizes -g ~/bwti/$geno/$geno.N-blocks.bed --avoid

fastaFromBed -name -fi ~/bwti/$geno/$geno.fa -bed random.bed -fo random.fa
perl -i -pe 's/(\w{50})/$1\n/g' random.fa

memeroot=../data/motifs/MEME
memedir=$memeroot/meme_out
fimodir=$memeroot/fimo_out
tomdir=$memeroot/tomtom_out
mkdir -p $memedir
mkdir $fimodir
mkdir $tomdir

#meme=/n/local/stage/meme/meme_4.8.1/src/meme.bin
#memep=/n/local/bin/meme/meme_p
meme=/n/local/stage/meme/meme_4.10.1/install/bin/meme
memep=/n/local/stage/meme/meme_4.10.1/install/bin/meme
fimo=/n/local/stage/meme/meme_4.10.1/install/bin/fimo
tomtom=/n/local/stage/meme/meme_4.10.1/install/bin/tomtom
if [ $cores -gt 1 ]; then 
    memecom="$memep -p $cores"
else
    memecom="$meme";
fi
nmotifs=20
peakfa=../data/macs2/$IP/$IP.macs2_peaks.top300.narrowPeak.fa
memecom="$memecom -dna -revcomp -mod zoops -nmotifs $nmotifs -minw 6 -maxw 15 -maxsize 10000000 -oc $memedir $peakfa"
echo $memecom
$memecom &

tomcom="$tomtom -oc $tomdir -min-overlap 4 $memedir/meme.html /n/data1/biobase/transfac/current/meme/transfac.meme";
echo $tomcom
$tomcom &

for i in {1..$nmotifs}
do
    fimodiri=$fimodir/$i
    mkdir $fimodiri
    fimofcom="$fimo --text --max-stored-scores 100000000 --motif $i $memedir/meme.html $peakfa > $fimodiri/fimo_fg.txt";
    fimobcom="$fimo --text --max-stored-scores 100000000 --motif $i $memedir/meme.html random.fa > $fimodiri/fimo_bg.txt";
    echo $fimofcom
    echo $fimobcom
    $fimofcom &
    $fimobcom &
done
