#!/usr/bin/env Rscript
library(methods)
library(GenomicRanges)
source("/n/projects/apa/R/apa_tools.R")

## Writes a correlation matrix to txt & png, for binned genome coverage vectors, found in RChipQC.RData objects
## Call "/home/apa/local/bin/ChIPBinCorr.R <outprefix> <sample1> <sample2> ... <sampleN>
##      where each <sample> has format "alias,/path/to/sample.RChipQC.RData"

## Setup

## Get command line
ca <- commandArgs(trailing=TRUE)
outprefix <- ca[1]
RDatas <- ca[2:length(ca)]   ## expecting "alias,/path/to/.RChipQC.RData"

## Parse sample strings
names(RDatas) <- sub(".RChipQC.RData$","",RDatas)
N <- length(RDatas)
RDatas2 <- matrix("", N, 2, TRUE, list(1:N,c("Alias","File")))
for (i in 1:N) {
    x <- unlist(strsplit(RDatas[i],"[=,]"))  # ideally, length-2 vector: alias, RData
    if (!file.exists(x[2])) stop(paste0("RChipQC RData file '",x[2],"' does not exist!  Stopping.\n"))
    RDatas2[i,] <- x
}

## Print sample table to screen
message("\nSample table:")
message(paste(colnames(RDatas2),collapse="\t"))
for (i in 1:N) message(paste(RDatas2[i,],collapse="\t"))

all.covg <- all.depth <- all.peaks <- all.stats <- new.list(RDatas2[,1])

for (i in 1:N) {
    true.i <- i
    message(paste(RDatas2[i,1],":",system("date",intern=T)))
    load(RDatas2[i,2])   ## will overwrite 'i'
    i <- true.i
    all.covg[[i]] <- attributes(bin.sums)$elementMetadata$sum
#    all.depth[[i]] <- rowSums(bpAtDepth)
#    all.peaks[[i]] <- allpeaks
#    all.stats[[i]] <- c(IP=ifelse(rdata[1]=="IP",1,0), read.len=readlen, genome=genome, aligns=aligns, unq.pos=upos, sequenced.bp=bpSeq, nonzero.pos=bpCov, unq.depths=L2, top1pct.pos=top1pct.pos, top1pct.seq=top1pct.seq)
#    suppressWarnings(rm(bin.sums,allpeaks,bpAtDepth,readlen,genome,aligns,upos,bpSeq,bpCov,L2,top1pct.pos,top1pct.seq))
}

all.covg <- do.call(cbind, all.covg)
#all.depth <- t(merge.table.list(all.depth))
#all.stats <- do.call(cbind,all.stats)

#cnas <- colnames(all.stats)

cmo <- corr.mat(log2(all.covg+1), plot=TRUE, main="Sample Correlations, Log2 Binned Genome Sums", imgdim=c(700,600), pmar=c(10,10), filename=paste0(outprefix,".corr.mat.ord.png"))
write.table2(cmo, paste0(outprefix,".corr.mat.ord.txt"), row.head="")

cmc <- corr.mat(log2(all.covg+1), plot=TRUE, reorder=TRUE, main="Sample Correlations, Log2 Binned Genome Sums", imgdim=c(700,600), pmar=c(10,10), filename=paste0(outprefix,".corr.mat.clust.png"))
write.table2(cmc, paste0(outprefix,".corr.mat.clust.txt"), row.head="")

message("ChIPBinCorr.R complete!\n")

## Stats from RChipQC.R 'stats.out' object
#        c("Bam file", bamfile2),
#        c("Read Length", readlen),
#        c("Min Height of Reportable Peak", minht),
#        c("Genome Pos", genome),
#        c("Nonzero Pos", bpCov),
#        c("Nonzero Pos %", round(100*bpCov/genome,2)),
#        c("Total Bp Aligned", bpSeq),
#        c("Pos in Top 1% of Depths", round(100*top1pct.pos,2)),
#        c("Bp in Top 1% Pos", top1pct.seq),
#        c("Bp in Top 1% Pos, as % Bp Aligned", round(100*top1pct.seq/bpSeq,2)),
#        c("Alignments", aligns),
#        c("Unique Alignment Starts", upos),
#        c("Unique Alignment %", round(100*upos/aligns,2)),
#        c("Total Islands", sum(sapply(islands,nrow))),
#        c("Islands >= Min Height", nrow(allpeaks)),
#        c("% Islands >= Min Height", round(100*nrow(allpeaks)/sum(sapply(islands,nrow)),2))


