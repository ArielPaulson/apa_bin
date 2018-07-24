#!/usr/bin/env Rscript

ca <- commandArgs(trailing=TRUE)
bamfile <- ca[1]
RData <- ca[2]

if (!file.exists(bamfile)) stop(paste0("BAM file '",bamfile,"' does not exist!\n"))
if (is.na(RData)) RData <- sub("\\.bam$",".RData",bamfile)

library(chipseq)
library(rtracklayer)

bam <- as(import(bamfile, "bam"),"GRanges")
covg <- coverage(bam)
save(covg, file=RData)
quit()

