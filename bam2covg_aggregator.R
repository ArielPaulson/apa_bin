#!/usr/bin/env Rscript

ca <- commandArgs(trailing=TRUE)
RData.out <- ca[1]
RData.in <- ca[2:length(ca)]
R <- length(RData.in)

fail <- FALSE
for (rd in RData.in) {
	if (!file.exists(rd)) {
	   message(paste0("Input RData object '",rd,"' does not exist!"))
	   fail <- TRUE
	}
}
if (fail) stop("Some inputs missing; halting.\n")


library(chipseq)
library(rtracklayer)

message(paste0("Initializing with '",RData.in[1],"'..."))
load(RData.in[1])  # 'covg'
cov.master <- covg

for (i in 2:R) {
	message(paste0("Adding '",RData.in[i],"'..."))
	load(RData.in[i])  # 'covg'
	L <- length(covg)
	for (j in 1:L) {
		w <- which(names(cov.master)==names(covg)[j])
		if (is.na(w)) {
		   stop(paste0("Input RData object ",i," ('",RData.i[i],"') has seqname '",names(covg)[j],"', but not the master object!\n"))
		} else {
		   cov.master[[w]] <- cov.master[[w]] + covg[[j]]
		}
	}
}

save(cov.master, file=RData.out)
message("bam2covg_aggregator.R complete!\n")
quit()

