## a wrapper function that fetches the available KEGG data and
## corresponding InterPro IDs from biomaRt for the universe of
## entrezgene IDs

getKEGGdata <- function(universe=NULL, pathways=NULL,
                        ensemblMart="hsapiens_gene_ensembl")
{
    ## check if we are online
    op <- options(warn=-1)
    on.exit(options(op))
    if(class(try(readLines("http://www.bioconductor.org"),
             silent=TRUE)) == "try-error")
        stop("Active internet connection needed for this function")
    options(op)
    ## get all available human KEGG pathway annotations and
    ## corresponding entrezgene IDs
    if(!is.null(pathways))
        hKEGGids <- pathways
    else
        hKEGGids <- grep("^hsa", ls(KEGGPATHID2EXTID), value=TRUE)
    path2Genes <- mget(hKEGGids, KEGGPATHID2EXTID)
    hKEGGgenes <- union(universe, unique(unlist(path2Genes, use.names=FALSE)))
    hKEGGgenes <-  hKEGGgenes[!is.na(hKEGGgenes)]

    ## get corresponding interpro domains from biomaRt
    ensembl <- useMart("ensembl", dataset=ensemblMart)
    tmp <- getBM(attributes=c("entrezgene", "interpro"), filters="entrezgene",
                 values=hKEGGgenes, mart = ensembl, bmHeader=FALSE)
    gene2Domains <- split(tmp$interpro, tmp$entrezgene, drop=FALSE)
    missing <- setdiff(hKEGGgenes, names(gene2Domains))
    gene2Domains[missing] <- ""
    hKEGGdomains <- unique(unlist(gene2Domains))
    hKEGGdomains <- hKEGGdomains[!is.na(hKEGGdomains)]
    path2Domains <- lapply(path2Genes, function(x, gene2Domains)
                           unique(unlist(gene2Domains[x], use.names=FALSE)), 
                           gene2Domains)
    
    ## the lengths of pathway, gene and domain vectors
    dims <- c(pathway=length(hKEGGids), gene=length(hKEGGgenes),
              domain=length(hKEGGdomains))
    
    ## create an object of class 'ipDataSource'
    return(new("ipDataSource", genes=hKEGGgenes, pathways=hKEGGids,
               domains=hKEGGdomains, gene2Domains=gene2Domains,
               path2Domains=path2Domains, dims=dims, type="KEGG"))
}

