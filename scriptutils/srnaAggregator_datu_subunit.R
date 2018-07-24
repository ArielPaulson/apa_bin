#!/usr/bin/env Rscript

source("/home/apa/apa_tools.R")

ca <- commandArgs(trailing=TRUE)
i <- as.numeric(ca[1])  # branch of 'dat2' to aggregate
s <- as.numeric(ca[2])  # sample name to aggregate
d <- ca[3]              # subunit ops dir
if (!grepl("/$",d)) d <- paste0(d,"/")

message("Loading prep.RData...")
load(paste0(d,"prep.RData"))

sub.datu <- lapply(bpi, function(b) {
    message(paste(i,s,b,":",names(AUi)[i],names(sampi)[s],names(bpi)[b]))
    if (length(dat2[[i]][[s]][[b]])>0) {
        x <- dat2[[i]][[s]][[b]]
        y <- ablank2
        if (any(x[,2] %in% y[,2])) {
            xa <- aggregate(x[,4:6], list(x[,2]), sum)
            y[,4:6] <- xa[match(y[,2],xa[,1]),2:4]
        }
        y
    }
})

save(sub.datu, file=paste0(d,"sub.datu.",i,".",s,".RData"))
message("srnaAggregator_datu_subunit.R complete!")


