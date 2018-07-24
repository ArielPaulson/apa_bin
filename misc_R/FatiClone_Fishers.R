
data <- read.delim("/home/apa/FatiClone_Fishers_input.txt", sep="\t", header=F, fill=T)
colnames(data) <- c("\tk","cluster","DB.level","GO.acc","cluster.in","cluster.out","bkg.in","bkg.out","clust.pct","bkg.pct","clust.enrich","p.raw","p.adj","odds","conf.int.lower","conf.int.upper")
nrow(data)
data2 <- as.matrix(data[,5:10])

fishers <- function (vec) {
	fmat <- matrix(data=vec[1:4], nrow=2, ncol=2, byrow=T)	# vec = c( cluster in, cluster out, bkg in, bkg out )
	altern <- 'two.sided'
	x <- fisher.test(fmat, alternative=altern, conf.lev=0.95)
	y <- c(x[[1]], x[[3]][[1]], x[[2]][1], x[[2]][2])	# raw p value, odds ratio, confidence interval lower bound, confidence interval upper bound
}

z <- as.matrix(apply(data2, 1, fishers))
data[,c(12,14,15,16)] <- t(z)
clusters <- unique(data[,2])
for (i in 1:length(clusters)) {
	cv <- which(data[,2] == clusters[i])
	data[cv,13] <- p.adjust(data[cv,12], method="BH")	# must adjust WITHIN "clusters" when using mono list, bkg=opposite
}


write.table(data, file="/home/apa/FatiClone_Fishers_ALL_output.txt", sep="\t", quote=F)
