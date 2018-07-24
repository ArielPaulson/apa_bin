#!/usr/bin/env Rscript

library(SPIA)
inpath <- commandArgs(trailing=TRUE)[[1]]
outpath <- commandArgs(trailing=TRUE)[[2]]
organism <- commandArgs(trailing=TRUE)[[3]]

makeSPIAdata(kgml.path=inpath, out.path=outpath, organism=organism)
quit()

