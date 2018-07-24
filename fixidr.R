#!/usr/bin/env Rscript

source("/home/apa/local/bin/apa_tools/code/math/fishers.method.R")
source("/home/apa/local/bin/apa_tools/code/IO/stabilize.number.R")
source("/home/apa/local/bin/apa_tools/code/IO/write.table2.R")

idrdir <- commandArgs(TRUE)[1]
if (!grepl("/$",idrdir)) idrdir <- paste0(idrdir,"/")
idr1 <- system(paste0("ls ",idrdir,"IDR1/*.expanded.txt 2>/dev/null"),intern=TRUE)
idr2 <- system(paste0("ls ",idrdir,"IDR2/*.expanded.txt 2>/dev/null"),intern=TRUE)

for (file in c(idr1,idr2)) {
    
    ##  1-14: Chrom,Start,End,Name,Score,Strand,FoldChange,Pvalue,Qvalue,Summit,Width,Height,IDR.Plocal,IDR.Pglobal
    ## 15-27: Pool.Chrom,Pool.Start,Pool.End,Pool.Name,Pool.Score,Pool.Strand,Pool.FoldChange,Pool.Pvalue,Pool.Qvalue,Pool.Summit,Pool.Width,Pool.Height,Pool.Rank
    ## 28-39: RepA.Chrom,RepA.Start,RepA.End,RepA.Name,RepA.Score,RepA.Strand,RepA.FoldChange,RepA.Pvalue,RepA.Qvalue,RepA.Summit,RepA.Width,RepA.Height
    ## 40-51: RepB.Chrom,RepB.Start,RepB.End,RepB.Name,RepB.Score,RepB.Strand,RepB.FoldChange,RepB.Pvalue,RepB.Qvalue,RepB.Summit,RepB.Width,RepB.Height
    
    message(file)
    x <- read.delim(file, as.is=TRUE)
    x$Pvalue <- -log10(apply(10^-as.matrix(x[,c("RepA.Pvalue","RepB.Pvalue")]), 1, fishers.method, zero=FALSE))
    x$Qvalue <- -log10(apply(10^-as.matrix(x[,c("RepA.Qvalue","RepB.Qvalue")]), 1, fishers.method, zero=FALSE))
    for (i in c(6,20,33,45)) x[[i]] <- "."   # correct all strands
    write.table2(x, file, stabilize=TRUE)
}
