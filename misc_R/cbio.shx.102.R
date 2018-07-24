
source("~/uapa/R/apa_tools.R")

x <- read.delim("~/bioan/RongLi/shx/cbio.shx.102/results/GOanalysis/SGvsCT/SGvsCT_sig_GO.txt", as.is=T)
x[,2] <- "SGvsCT"
y <- read.delim("~/bioan/RongLi/shx/cbio.shx.102/results/GOanalysis/STvsCT/STvsCT_sig_GO.txt", as.is=T)
y[,2] <- "STvsCT"
z <- read.delim("~/bioan/RongLi/shx/cbio.shx.102/results/GOanalysis/overlaps/sig_overlaps_GO.txt", as.is=T)
z[,2] <- "overlap"


w <- rbind(x,y,z)
w[,1] <- 1

go.hm <- GO.heatmap(w[,c(1,2,6,13)], ret.data=T)
