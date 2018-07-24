#!/usr/bin/env Rscript

## Last run on 20161118

system("wget -O unicode_names_raw.txt http://unicode.org/Public/UNIDATA/NamesList.txt")

x <- scan("unicode_names_raw.txt", what="", sep="\n"); length(x)
y <- x[grepl("^[0-9A-Z]",x)]; length(y)
z <- as.data.frame(do.call(rbind, strsplit(y,"\t"))); dim(z)
w <- cbind(strtoi(paste0("0x",z[,1])),z)
colnames(w) <- c("Dec","Hex","Name")
write.table(w, "unicode_names.txt", sep="\t", quote=FALSE, row.names=FALSE)


