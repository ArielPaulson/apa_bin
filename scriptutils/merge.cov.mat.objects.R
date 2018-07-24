
## be in the cov.mats directory you want to merge into

other.dir.full <- "../cov.mats.patch"   # directory to merge from
other.dir <- sub(".*/","",other.dir.full)
this.dir <- sub(".*/","",getwd())

## Load other data
load(paste0(other.dir.full,"/cov.mat.beds.RData"));  b=cov.mat.beds
load(paste0(other.dir.full,"/cov.mat.paths.RData")); p=cov.mat.paths

## Load this data
load("cov.mat.beds.RData")
load("cov.mat.paths.RData")

this.len <- length(cov.mat.beds)   # either beds or paths will do
other.len <- length(b)
this.len; other.len

names(cov.mat.beds);  names(b)
names(cov.mat.paths); names(p)






## DIFFERENT-BEDS WORKFLOW

cov.mat.beds  <- c(cov.mat.beds, b)
cov.mat.paths <- c(cov.mat.paths, p)

names(cov.mat.beds)
names(cov.mat.paths)

other.pos <- c(1:other.len)+this.len

for (o in other.pos) {
    k.list <- grep("^k",names(cov.mat.paths[[o]]))
    for (k in k.list) {   # k5, k2, etc
        for (i in 1:length(cov.mat.paths[[o]][[k]])) {  # IP, input
            cov.mat.paths[[o]][[k]][[i]] <- sub(other.dir,this.dir,cov.mat.paths[[o]][[k]][[i]])
        }
    }
}

save(cov.mat.beds, file="cov.mat.beds.RData")
save(cov.mat.paths, file="cov.mat.paths.RData")






## SAME-BEDS WORKFLOW
## **** Assuming k-set is equal ****

for (n in names(cov.mat.paths)) {
    t.n <- which(names(cov.mat.paths)==n)
    o.n <- which(names(p)==n)
    message(t.n,",",o.n)
    k.list <- grep("^k",names(cov.mat.paths[[t.n]]))
    for (k in k.list) {   # k5, k2, etc
        for (i in 1:length(cov.mat.paths[[t.n]][[k]])) {  # IP, input
            message(t.n,",",o.n,":",k,",",i,":",length(p[[o.n]][[k]]))
            if (length(p[[o.n]][[k]])>=i) {
                if (length(p[[o.n]][[k]][[i]])>0) {
                    add <- sub(other.dir,this.dir,p[[o.n]][[k]][[i]])
                    add <- add[!(names(add) %in% names(cov.mat.paths[[t.n]][[k]][[i]]))]
                    message(t.n,",",o.n,":",k,",",i,":",length(p[[o.n]][[k]]),",",length(add))
                    cov.mat.paths[[t.n]][[k]][[i]] <- c(cov.mat.paths[[t.n]][[k]][[i]],add)
                }
            }
        }
    }
}

save(cov.mat.paths, file="cov.mat.paths.RData")

