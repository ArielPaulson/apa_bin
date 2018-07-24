#!/usr/bin/env Rscript

source("/home/apa/apa_tools.R")
source("/home/apa/local/bin/scriptutils/snpOrthology_functions.R")


RData <- commandArgs(TRUE)[1]   # 'finalize.RData' output from snpOrthology
load(RData)

ortho <- aligns <- exons <- vars <- pdb.report <- new.list(genes[,1])

lost1 <- lost2 <- rep(FALSE,nrow(genes))


for (i in 1:G) {
    message(i,"/",G," ",genes[i,1])
    if (genes[i,1] %in% no.ortho) next
    rdatai <- paste0(resultsdir,genes[i,1],".RData")
    if (!file.exists(rdatai)) {
        message("WARNING: expected RData file '",rdatai,"' does not exist!")
        lost1[i] <- TRUE
    } else {
        load(rdatai)  # 'x'
        if (length(x)>0) {
            vars[[i]] <- x$vars2   # updated
            ortho[[i]] <- x$ortho  # updated
            exons[[i]] <- x$exons
            aligns[[i]] <- x$aligns
            if (do.PDB & genes[i,1] %in% pdbs$Transcript) pdb.report[[i]] <- x$report
        } else {
            message("WARNING: expected RData file '",rdatai,"' has no contents!")
            lost2[i] <- TRUE
        }
    }
}


if (any(lost1)) message(sum(lost1)," transcripts did not complete!  This may or may not be OK.")
if (any(lost2)) message(sum(lost2)," transcripts do not have results!  This may or may not be OK.")
lost3 <- rbind2( cbind(1,genes[lost1,1]), cbind(2,genes[lost2,1]) )
write.table2(lost3, paste0(outprefix,"lost.txt"), col.names=("Loss.Type","Transcript.ID"))


## Prepare and output PDB reports
if (do.PDB) {
    message("Compiling PDB report...")
    pdb.report2 <- pdb.report
    for (i in 1:length(pdb.report)) {
        for (p in 1:length(pdb.report[[i]])) {
            x <- do.call(rbind, lapply(1:length(pdb.report[[i]][[p]]), function(j) suppressWarnings(cbind(Variant=names(pdb.report[[i]][[p]])[j], pdb.report[[i]][[p]][[j]])) ))
            pdb.report2[[i]][[p]] <- x[order(x[,7]),]
            pdb.report2[[i]][[p]][,7] <- round(pdb.report2[[i]][[p]][,7],3)
        }
        pdb.report2[[i]] <- do.call(rbind, lapply(1:length(pdb.report2[[i]]), function(p) suppressWarnings(cbind(PDB=names(pdb.report2[[i]])[p], pdb.report2[[i]][[p]])) ))
    }
    pdb.report2 <- rownameless(do.call(rbind, lapply(1:length(pdb.report2), function(i) suppressWarnings(cbind(Gene=genes[genes[,1]==names(pdb.report2)[i],2], pdb.report2[[i]])) )))
    pdb.report2 <- pdb.report2[,c(1,3,2,4:ncol(pdb.report2))]
    colnames(pdb.report2) <- sub("\\.1$",".NBR",colnames(pdb.report2))
    write.table(pdb.report2, paste0(outprefix,"PDB_report.txt"), sep="\t", na="", quote=FALSE, row.names=FALSE)
}


## Write 'vars' tables and exit
message("Compiling variation table...")
llv <- listLengths(vars,nrow)
maxdf <- vars[[which.max(listLengths(vars,ncol))]][0,]
vars3 <- rownameless(maxdf[1:sum(llv),])
row <- 0
for (i in which(llv>0)) {
    x <- vars[[i]]
    nrx <- nrow(x)
    vars3[(row+1):(row+nrx),match(colnames(x),colnames(vars3))] <- x
    row <- row + nrx
}

write.table(vars3, paste0(outprefix,"variation_table.txt"), sep="\t", na="", quote=FALSE, row.names=FALSE)
save.image(RData)
quit()

