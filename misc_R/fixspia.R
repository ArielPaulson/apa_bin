#!/usr/bin/env Rscript

rdata <- commandArgs(trailing=TRUE)[1]
rdata2 <- sub(".RData","-no-Olfactory.RData",rdata)
load(rdata)
is.olfac <- names(path.info) %in% "04740"
if (any(is.olfac)) {
    path.info <- path.info[!is.olfac]
    save(path.info, file=rdata2)
}
