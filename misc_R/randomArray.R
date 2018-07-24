
x <- as.numeric(commandArgs(trailingOnly=TRUE))
M <- x[4]-x[3]+1

y <- roundc(runif(x[2], x[3], x[4]), 0)

if (!as.logical(x[5])) {
   if (M >= x[2]) {
      dups <- duplicated(y)
      duped <- any(dups)
      while (duped) {
      	  y[dups] <- round(runif(sum(dups), x[3], x[4]), 0)
      	  dups <- duplicated(y)
      	  duped <- any(dups)
      }
      write(y, file=paste(c("randomArray",x,"tmp"), collapse="."), ncol=1)
   } else {
      message(paste("Cannot sample",x[2],"values from",x[3],"-",x[4],"without replacement!"))
   }
} else {
   write(y, file=paste(c("randomArray",x,"tmp"), collapse="."), ncol=1)
}

