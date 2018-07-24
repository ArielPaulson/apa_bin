
afiles <- system("ls *hrs*.htseq.counts.txt", intern=T)
names(afiles) <- sub(".htseq.*","",afiles)
alist <- lapply(afiles, read.delim, as.is=T, header=F)
adat <- as.matrix(do.call(cbind, lapply(alist, "[[", 2)))[,c(7:10,1:6)]
rownames(adat) <- alist[[1]][,1]
adat <- adat[!grepl("_",rownames(adat)),]
rm(alist)

bfiles <- system("ls *.htseq.counts.txt | grep -v hrs", intern=T)
names(bfiles) <- sub(".htseq.*","",bfiles)
blist <- lapply(bfiles, read.delim, as.is=T, header=F)
bdat <- as.matrix(do.call(cbind, lapply(blist, "[[", 2)))[,c(2,5,1,4,3,6)]
rownames(bdat) <- blist[[1]][,1]
bdat <- bdat[!grepl("_",rownames(bdat)),]
rm(blist)

peaks <- read.delim("../../data/Hoxa1_old_vs_new/FINAL_REFILTERED_HOXA1_PEAKLIST.bed", as.is=T, header=F)
nns <- read.delim("../../data/Hoxa1_old_vs_new/FINAL_REFILTERED_HOXA1_PEAKLIST.nn.txt", as.is=T)





library(edgeR)
gdat <- read.delim("~/bwti/mm10/Ens_80/mm10.Ens_80.genedata.txt",as.is=T)
gdat <- gdat[order(gdat[,1]),]
bio <- gdat$SimpleBiotype %in% qw(protein_coding,lncRNA,pseudogene)

all(rownames(adat)==rownames(bdat))
all(rownames(adat)==gdat[,1])

e2m <- read.delim("../ens80_mgi.txt",as.is=T)
e2m2 <- breakout(lapply(split(e2m[,2],e2m[,1]), function(x) paste(suniq(x),collapse=";")),rev=T)

group <- list(A=factor2(sub("_RA.*","",colnames(adat))), B=factor2(sub("_[12]","",colnames(bdat))))
lapply(group,as.numeric)
grpi <- c(A=1,B=2)
G <- length(group)
pairs <- list(A=list(`6-4`=2:1,`12-6`=3:2,`24-12`=4:3,`36-24`=5:4), B=list(`Dox-Un`=2:1,`DoxRA-Dox`=3:2,`DoxRA-Un`=c(3,1)))  # comparisons to be made
ecounts <- list(A=adat[bio&apply(adat>0,1,any),], B=bdat[bio&apply(bdat>0,1,any),])
egdat <- list(A=gdat[match(rownames(ecounts$A),gdat[,1]),], B=gdat[match(rownames(ecounts$B),gdat[,1]),])
sapply(ecounts,nrow); sapply(egdat,nrow)
dge <- lapply(grpi, function(i) estimateTagwiseDisp( estimateCommonDisp( calcNormFactors( DGEList(counts=ecounts[[i]], group=group[[i]]) ) ) ) )

et <- top <- sig <- to.GO <- pairs  # have i, j components
norm.count <- all.rpkm <- group   # have i component only
all.xls <- sig.xls <- hmx <- hmz <- hcx <- hcz <- c(pairs[[1]],pairs[[2]])  # have j component only
to.GO.ens <- data.frame(GENE="",GROUP="")[0,]

for (i in 1:G) {
    ## DO THESE SEPARATELY FIRST
    norm.count[[i]] <- list(
        mean=aggregate.cols(dge[[i]]$pseudo.counts, list(group[[i]]), mean),
        sd=aggregate.cols(dge[[i]]$pseudo.counts, list(group[[i]]), sd),
        cv=0,
        rpkm=0
    )
    colnames(norm.count[[i]]$mean) <- paste0(colnames(norm.count[[i]]$mean), ".CpmAvg")
    colnames(norm.count[[i]]$sd) <- paste0(colnames(norm.count[[i]]$sd), ".CpmSd")
    norm.count[[i]]$cv <- norm.count[[i]]$sd/norm.count[[i]]$mean
    norm.count[[i]]$rpkm <- 1E3*norm.count[[i]]$mean/egdat[[i]]$Uxon_Len   # 1E3 NOT 1E9 !!!
    for (j in 1:length(norm.count)) norm.count[[i]][[j]] <- norm.count[[i]][[j]][match(egdat[[i]][,1],rownames(norm.count[[i]][[j]])),]
    all.rpkm[[i]] <- 1E3*dge[[i]]$pseudo.counts/egdat[[i]]$Uxon_Len   # 1E3 NOT 1E9 !!!
}

k <- 0
for (i in 1:G) {
    for (j in 1:length(pairs[[i]])) {
        k <- k+1
        
        FDR <- ifelse(i==1&j==1,0.05,0.001)
        IM(names(pairs)[i],names(pairs[[i]])[j])
        et[[i]][[j]] <- exactTest(dge[[i]], pair=pairs[[i]][[j]])
        et[[i]][[j]]$table <- cbind(et[[i]][[j]]$table, decideTestsDGE=decideTestsDGE(et[[i]][[j]], p=FDR, adjust="BH"))
        top[[i]][[j]] <- topTags(et[[i]][[j]], nrow(ecounts[[i]]))
        top[[i]][[j]]$table <- top[[i]][[j]]$table[,c(1:3,5,4)]
        colnames(top[[i]][[j]]$table)[5] <- "DE"
        sig[[i]][[j]] <- top[[i]][[j]]$table[top[[i]][[j]]$table$DE!=0,]
        IM(names(pairs)[i],names(pairs[[i]])[j],sum(sig[[i]][[j]]$DE==1),sum(sig[[i]][[j]]$DE==-1))
        
        to.GO[[i]][[j]] <- cbind(GENE=rownames(sig[[i]][[j]]),GROUP=paste0(sub("1","",sign(sig[[i]][[j]]$DE)),names(pairs[[i]][j])))
        to.GO.ens <- rbind(to.GO.ens, to.GO[[i]][[j]])
        
        tt <- top[[i]][[j]]$table
        eg <- egdat[[i]][match(rownames(tt),egdat[[i]][,1]),]
        rpkm <- all.rpkm[[i]]
        colnames(rpkm) <- paste0(colnames(rpkm),".RPKM")
        all.xls[[k]] <- cbind(eg[,1:2], tt, norm.count[[i]]$mean, norm.count[[i]]$sd, rpkm, eg[,3:ncol(eg)])
        sig.xls[[k]] <- all.xls[[i]][all.xls[[i]]$DE!=0,]
        
        hmx[[k]] <- as.matrix(log2(all.rpkm[[i]][rownames(all.rpkm[[i]]) %in% rownames(sig[[i]][[j]]),]+1))
        colnames(hmx[[k]]) <- sub(".RPKM","",colnames(hmx[[k]]))
        hmz[[k]] <- row.norm(hmx[[k]],divSD=T)
        hcx[[k]] <- list(
            row=reorder.hclust2(hclust(dist(hmx[[k]]),"average"), hmx[[k]], mean)$order,
            col=hclust(dist(t(hmx[[k]])),"average")$order
        )
        hcz[[k]] <- list(
            row=hclust(distance(hmz[[k]]),"average")$order,
            col=hclust(distance(t(hmz[[k]])),"average")$order
        )
        
        png(paste0(names(pairs[[i]])[j],".hmx.png"), 1000, max(c(500,nrow(hmx[[k]]))))
        mipu(hmx[[k]][hcx[[k]]$row,], pal="Reds", cex=1.2, pmar=c(10,1), main=paste(names(pairs[[i]])[j],"DE Gene Log2(RPKMs+1)"))
        dev.off()
        png(paste0(names(pairs[[i]])[j],".hmz.png"), 1000, max(c(500,nrow(hmx[[k]]))))
        mipu(hmz[[k]][hcz[[k]]$row,], pal="RWB", col.center=0, col.limits=c(-4,4), cex=1.2, pmar=c(10,1), main=paste(names(pairs[[i]])[j],"DE Gene Z-Scores"))
        dev.off()
    }
}

to.GO.mgi <- to.GO.ens
to.GO.mgi[,1] <- e2m2[match(to.GO.ens[,1],e2m2[,1]),2]
write.table2(to.GO.ens, "to.GO.ens.txt")
write.table2(to.GO.mgi, "to.GO.mgi.txt")

#WriteXLS2(all.xls, "all.edgeR.genes.xls", BoldHeaderRow=T, row.names=F, FreezeRow=1, FreezeCol=2)
WriteXLS2(sig.xls, "sig.edgeR.genes.xls", BoldHeaderRow=T, row.names=F, FreezeRow=1, FreezeCol=2)




save.image("Feb18.RData")



