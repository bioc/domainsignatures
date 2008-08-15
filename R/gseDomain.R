## The main function to run geneset enrichment based on Interpro
## domain signatures. It takes an ipDataSource object and the gene list
## as entrezgene IDs.

gseDomain <- function(dataSource, geneset, n=10000, verbose=TRUE,
                      samples=FALSE){
  ## validate arguments
  missing <- setdiff(geneset, dataSource@genes)
  if(length(missing)>0)
    stop("There are genes missing from the universe:\n",
         paste("   ", missing,collapse="\n"))

  ## start computation
  res <- compSimilarities(geneset, dataSource, n=n, verbose=verbose,
                          samples=samples)
  return(res)
}
