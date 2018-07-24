#!/usr/bin/env Rscript

source("/home/apa/local/bin/apa_tools/code/IO/read.macs2.xls.R")
source("/home/apa/local/bin/apa_tools/code/IO/read.narrowPeak.R")
source("/home/apa/local/bin/apa_tools/code/IO/write.bed.R")
source("/home/apa/local/bin/apa_tools/code/IO/stabilize.number.R")

outpref <- sub("macs2_peaks.*","",commandArgs(TRUE)[1])   # input is either a narrowPeak file, or the peakset-rep prefix
if (!grepl("\\.$",outpref)) outpref <- paste0(outpref,".")

xls <- paste0(outpref,"macs2_peaks.xls.gz")
xls1 <- sub(".gz$","",xls)
system(paste0("zcat ",xls," | grep \"^#\" > ",xls1))
system(paste0("zcat ",xls," | grep -v \"^#\" | grep -P \"\\S\" | sed 's/ //g' >> ",xls1))
system(paste0("gzip -f ",xls1))
xl <- read.macs2.xls(xls)

np1 <- paste0(outpref,"macs2_peaks.narrowPeak.gz")
np <- read.narrowPeak(np1)
xl <- xl[match(np[,4],xl[,10]),]         # xl in np1 order
np[,5] <- xl$pileup                      # replace score with peak height
write.bed(np, np1, format="narrowPeak")  # overwrite np1 file

np2 <- paste0(outpref,"macs2_peaks.sort100k.narrowPeak")
if (file.exists(np2)) {
    nps <- read.narrowPeak(np2)
    xl <- xl[match(nps[,4],xl[,10]),]         # xl in nps order (do this last; nps is a subset of np1)
    nps[,5] <- xl$pileup                      # replace score with peak height
    write.bed(nps, np2, format="narrowPeak")  # overwrite np2 file
}
