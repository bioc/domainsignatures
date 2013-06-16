## a function that takes the mapping between genes and arbitrary
## pathways as input (a list), fetches the necessary InterPro domains from
## biomaRt and creates an 'ipDataSource' object

dataSource <- function(mapping, type="generic"){
  if(is.null(names(mapping)))
    stop("'mapping' must be named list of pathway mappings to genes\n")
  pathways <- as.character(unique(unlist(mapping)))
  genes <- as.character(names(mapping))

  ## get corresponding interpro domains from biomaRt
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  tmp <- getBM(attributes=c("entrezgene", "interpro"), filters="entrezgene",
               values=genes, mart = ensembl, bmHeader=FALSE)
  gene2Domains <- split(tmp$interpro, tmp$entrezgene, drop=FALSE)
  missing <- setdiff(genes, names(gene2Domains))
  gene2Domains[missing] <- ""
  domains <- unique(unlist(gene2Domains))
  domains <- domains[!is.na(domains)]
  tmp2 <- unlist(mapping, use.names = FALSE)
  path2Genes <- split(rep(genes, listLen(mapping)), tmp2)
  path2Domains <- lapply(path2Genes, function(x, gene2Domains) unique(unlist(gene2Domains[x])), 
                         gene2Domains)

  ## the lengths of pathway, gene and domain vectors
  dims <- c(pathway=length(pathways), gene=length(genes),
            domain=length(domains))

  
  ## create an object of class 'ipDataSource'
  return(new("ipDataSource", genes=genes, pathways=pathways,
             domains=domains, gene2Domains=gene2Domains,
             path2Domains=path2Domains, dims=dims, type=type))
}
