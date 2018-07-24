x=as.numeric(system("cut -f3 -d' ' ../data/alignToRNA/H12Arep1/H12Arep1.trimmomatic.log", intern=T))
dhist(sample(x,100000))

y=table(x)
z=cbind(as.numeric(names(y)),c(y))
zf=zerofy(z[match(1:max(z[,1]),z[,1]),])
zf[,1]=1:max(z[,1])
lineplot(zf[,2])


x=read.delim("../data/alignToRNA/H12Arep1/H12Arep1.trim.histo.txt")
z=table(x[,1:2])
rnz <- as.numeric(rownames(z))
zf=zerofy(z[match(1:max(rnz),rnz),])
zf[,1]=1:max(rnz)
dhist(log10(zf))
