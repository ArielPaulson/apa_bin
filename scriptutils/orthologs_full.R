



########## NOT READY YET ##########




## setup
library(biomaRt)
sequence <- FALSE  # fetch ortholog sequences?
gene <- "ENSDARG00000034195"  # Danio
orgs <- c("Homo sapiens","Gallus gallus","Drosophila melanogaster","Danio rerio")  # SOURCE ORGANISM MUST BE LAST
orgs.ens <- sapply(orgs, function(x) paste0(tolower(substr(x,1,1)),unlist(strsplit(x," "))[2]) ); orgs.ens   # ensemblized organism names
N <- length(orgs)
M <- N-1  # everything but source org


## establish mart connections for each organism mart
marts <- as.list(orgs.ens); names(marts) <- orgs.ens
for (i in 1:N) marts[[i]] <- useMart("ensembl",paste0(orgs.ens[i],"_gene_ensembl"))  # create mart connections
la <- lapply(marts, listAttributes)  # available attributes for each mart


## find pertinent *ology fields in each mart
source2ortho.fields <- lapply(1:N, function(i) {
    ## attributes from source-org mart for all ortholog orgs, including self
    ## "log_" targeting "(ortho|para|homo)log_" attribute names
    grep("log_",grep(orgs.ens[i],la[[N]][,1],value=T),value=T)
})
names(source2ortho.fields) <- orgs.ens
ortho2source.fields <- lapply(1:M, function(i) {  # attributes from ortho-org marts for source org
    ## attributes from ortho-org marts for source org
    ## "log_" targeting "(ortho|para|homo)log_" attribute names
    ## strip organism string (e.g. "hsapiens_") from beginning of attribute name, to see which fields are common for all ortho orgs
    grep("log_",grep(orgs.ens[N],la[[i]][,1],value=T),value=T)
})
names(ortho2source.fields) <- orgs.ens[1:M]
source2ortho.fields
ortho2source.fields


## exhaustively detect pertinent orthologs, since many ortholog annotations are not bidirectional
source2ortho.genes <- lapply(1:N, function(i) {
    ## attributes from source-org mart for all ortholog orgs, including self (paralog test)
    ## "log_" targeting "(ortho|para|homo)log_" attribute names
    getBM(attributes=c(source2ortho.fields[[i]][1],"ensembl_gene_id"), filters="ensembl_gene_id", values=gene, mart=marts[[N]])
})
names(source2ortho.genes) <- orgs.ens
ortho2source.genes <- lapply(1:M, function(i) {  # attributes from ortho-org marts for source org
    ## attributes from ortho-org marts for source org
    ## "log_" targeting "(ortho|para|homo)log_" attribute names
    ## strip organism string (e.g. "hsapiens_") from beginning of attribute name, to see which fields are common for all ortho orgs
    message(orgs[i])
    x <- getBM(attributes="ensembl_gene_id", filters=paste0("with_homolog_",substr(orgs.ens[N],1,4)), values=TRUE, mart=marts[[i]])  # all genes with any source-org orthology
    y <- getBM(attributes=c("ensembl_gene_id",ortho2source.fields[[i]][1]), filters="ensembl_gene_id", values=x, mart=marts[[i]])  # gene-ortholog pairs
    y[y[,2]==gene,]  # gene-ortholog pairs where ortholog is the target gene
})
names(ortho2source.genes) <- orgs.ens[1:M]
source2ortho.genes
ortho2source.genes
#sapply(source2ortho.genes,nrow); t(sapply(source2ortho.genes, function(x) apply(x,2,function(y) length(unique(y)))))
#sapply(ortho2source.genes,nrow); t(sapply(ortho2source.genes, function(x) apply(x,2,function(y) length(unique(y)))))


## create master matrix with all ortholog annotations
xlogify <- function(x) unique(gsub("[a-z]+log(y?)_","xlog\\1_",sub("^[^_]+_","",x)))  # convert all (homo|ortho|para)log labels to "xlog", for unified accounting
cn1 <- c("Mart.Org","Ortho.Org")  # initial colnames for matrix
cn2 <- c("ensembl_gene_id",xlogify(unlist(c(source2ortho.fields,ortho2source.fields))))  # further colnames for matrix -- mart *log field types, with org prefixes removed.  

cn3 <- c("xlog_xlogy_type","xlog_subtype","xlog_xlogy_confidence","xlog_perc_id","xlog_perc_id_r1","xlog_dn xlog_ds")  # data unique to xlogs -- can get the rest from gene-data calls

source2ortho.rows <- lapply(1:N, function(i) {
    ## like source2ortho.genes, but returning all orthology fields
    x <- getBM(attributes=c("ensembl_gene_id",source2ortho.fields[[i]]), filters="ensembl_gene_id", values=gene, mart=marts[[N]])
    ncx <- ncol(x)
    colnames(x)[2:ncx] <- xlogify(colnames(x)[2:ncx])
    x <- x[,match(cn2,colnames(x))]  # match to output matrix column order
    colnames(x)[1:2] <- c("Target.Gene","Ortho.Gene")
    data.frame(Mart.Org=orgs[N],Ortho.Org=orgs[i],x)
})
names(source2ortho.rows) <- orgs.ens
ortho2source.rows <- lapply(1:M, function(i) {
    ## like ortho2source.genes, but returning all orthology fields
    x <- getBM(attributes=c("ensembl_gene_id",ortho2source.fields[[i]]), filters="ensembl_gene_id", values=ortho2source.genes[[i]][,1], mart=marts[[i]])
    ncx <- ncol(x)
    colnames(x)[2:ncx] <- xlogify(colnames(x)[2:ncx])
    x <- x[,match(cn2,colnames(x))]  # match to output matrix column order
    colnames(x)[1:2] <- c("Target.Gene","Ortho.Gene")
    data.frame(Mart.Org=orgs[i],Ortho.Org=orgs[N],x)  #[,c(1:2,4:3,5:(ncol(x)+2))]
})
names(ortho2source.rows) <- orgs.ens[1:M]
source2ortho.rows
ortho2source.rows

master <- do.call(rbind, c(source2ortho.rows,ortho2source.rows)); master


## add sequence, if requested
if (sequence) {
    
    
    
}

