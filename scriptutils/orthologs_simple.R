
## simple ortholog retrieval, given an ensembl gene ID and a list of target organisms

## setup
library(biomaRt)
sequence <- FALSE  # fetch ortholog sequences?
martname <- "ensembl"  # **** BASE ENSEMBL ONLY: if using EnsemblGenomes, use listMarts() to see alternate marts.
#martname <- "fungi_mart_19"  # example EnsemblGenomes mart name for fungi.ensembl.org.  Use pattern (fungi|metazoa|plants|protists)_mart_(version#).
gene <- "ENSDARG00000034195"  # Danio
orgs <- c("Homo sapiens","Gallus gallus","Drosophila melanogaster","Danio rerio")  # SOURCE ORGANISM MUST BE LAST
orgs.ens <- sapply(orgs, function(x) paste0(tolower(substr(x,1,1)),unlist(strsplit(x," "))[2]) ); orgs.ens   # ensemblized organism names
N <- length(orgs)
M <- N-1  # everything but source org


## establish mart connections for each organism mart
marts <- as.list(orgs.ens); names(marts) <- orgs.ens
for (i in 1:N) marts[[i]] <- useMart(martname,paste0(orgs.ens[i],"_gene_ensembl"))  # create mart connections
la <- lapply(marts, listAttributes)  # available attributes for each mart


## find pertinent orthology fields in each mart
source2ortho.fields <- lapply(1:M, function(i) {
    ## attributes from source-org mart for all ortholog orgs
    ## "log_" targeting "(ortho|para|homo)log_" attribute names
    grep("log_",grep(orgs.ens[i],la[[N]][,1],value=T),value=T)
})
names(source2ortho.fields) <- orgs.ens[1:M]


## get ortholog genes
source2ortho.genes <- lapply(1:M, function(i) {
    ## attributes from source-org mart for all ortholog orgs, including self (paralog test)
    ## "log_" targeting "(ortho|para|homo)log_" attribute names
    getBM(attributes=c(source2ortho.fields[[i]][1],"ensembl_gene_id"), filters="ensembl_gene_id", values=gene, mart=marts[[N]])
})
names(source2ortho.genes) <- orgs.ens[1:M]


## get ortholog data
source2ortho.data <- lapply(1:M, function(i) {
    ## like source2ortho.genes, but returning all orthology fields
    getBM(attributes=c("ensembl_gene_id",source2ortho.fields[[i]]), filters="ensembl_gene_id", values=gene, mart=marts[[N]])
})
names(source2ortho.data) <- orgs.ens[1:M]


## get ortholog sequence, if desired (adds new column)
if (sequence) {
    for (i in 1:M) {
        ## like source2ortho.genes, but returning only sequence
        col <- grep("log_ensembl_peptide",colnames(source2ortho.data[[i]]))[1]  # find the peptide ID column ([1] just in case?)
        source2ortho.data[[i]] <- data.frame(
            source2ortho.data[[i]],
            peptide=getSequence(id=source2ortho.data[[i]][,col], type="ensembl_peptide_id", seqType="peptide", mart=marts[[i]])$peptide
        )
    }
}

## extra steps which convert the list to a single data.frame
xlogify <- function(x) unique(gsub("[a-z]+log(y?)_","xlog\\1_",sub("^[^_]+_","",x)))  # convert all (homo|ortho|para)log labels to "xlog", for a unified accounting
master <- do.call(rbind, lapply(source2ortho.data, function(x){ y=x; colnames(y)=xlogify(colnames(y)); y }))

write.table(master, paste0(gene,".orthologs.txt"), sep="\t", quote=F, row.names=F)
save.image(paste0(gene,".orthologs.RData"))

