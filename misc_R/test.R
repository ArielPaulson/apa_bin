#!/bin/env Rscript

x <- matrix(rnorm(10000), nrow=1000)
plot(1:10, x[1,], type="l")
system.time({ for (i in 2:1000) lines(1:10, x[i,]) })


