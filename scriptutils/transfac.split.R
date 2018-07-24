#!/usr/bin/env Rscript
source("~/apa_tools.R")

ver <- commandArgs(TRUE)[1]

path <- paste0("/n/data1/biobase/transfac/",ver,"/meme/")
meme <- read.minimal.meme(paste0(path,"transfac.",ver,".meme"))

clades <- c("vertebrate","plant","fungal","insect","nematode")
prefs <- c("V","P","I","F","N")
matched <- rep(FALSE, length(meme$MOTIFS))

for (i in 1:5) {
    in.clade <- grep(paste0(" ",prefs[i],"\\$"), names(meme$MOTIFS))
    message(clades[i]," ",length(in.clade))
    m <- meme
    m$MOTIFS <- m$MOTIFS[in.clade]
    matched[in.clade] <- TRUE
    write.minimal.meme(m, paste0(path,"transfac.",ver,".",clades[i],".meme"))
    system(paste0("ln -sf ",path,"transfac.",ver,".",clades[i],".meme ",path,"transfac.",clades[i],".meme"))
}

message(sum(matched),"/",length(matched))

if (any(!matched)) {
    in.misc <- which(!matched)
    m <- meme
    m$MOTIFS <- m$MOTIFS[in.misc]
    matched[in.misc] <- TRUE
    write.minimal.meme(m, paste0(path,"transfac.",ver,".misc.meme"))
    system(paste0("ln -sf ",path,"transfac.",ver,".misc.meme ",path,"transfac.misc.meme"))
}

