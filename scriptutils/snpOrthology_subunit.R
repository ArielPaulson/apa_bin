#!/usr/bin/env Rscript

load(commandArgs(TRUE)[1])   # '/datadir/<transcript_ID>.RData' output from snpOrthology

source("/home/apa/apa_tools.R")
source("/home/apa/local/bin/scriptutils/snpOrthology_functions.R")

system(paste("mkdir -p",this.tmp.dir))
system(paste0("touch ",this.tmp.prefix,"dead"))
system(paste("rm -f",this.result.rdata))
x <- align.orthologs()
if (do.PDB) if (length(this.trans.pdbs)>0) x <- process.pdb(x)
save(x, file=this.result.rdata)
system(paste0("rm -f ",this.tmp.prefix,"dead"))
quit()
