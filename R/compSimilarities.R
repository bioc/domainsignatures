## For all pathways in an 'ipDataSource' object, compute chisquare-based
## similarity of a gene list to each pathway and provide subsampling
## distributions and p-values. Output of the function is a
## list with items:
##    dist: a list containing the subsampling distributions
##    similarity: a list of distance measures of the gene
##                  list to each pathway
##    pvalue: a named vector of p-values indicating similarity of
##            gene list to each pathway
compSimilarities <- function(geneList, dataSource, n=10000, verbose=TRUE,
                             samples=FALSE){
  geneList <- geneList[geneList %in% dataSource@genes]
  np <- length(dataSource@pathways)
  dists <- setDists <- list()
  ## setup progress report
  if(verbose && capabilities("tcltk")){
    require(prada)
    mess <- paste("resampling gene list\nwith ", n ,"samples")
    sub <- paste("(pathway 1 of ", np, ")", sep="")
    progress(message=mess, sub=sub)
    on.exit(killProgress())
  }
  ## all domains of the geneList
  listDoms <- unique(unlist(dataSource@gene2Domains[geneList], 
                     use.names=FALSE))
  ## iterate over all pathways
  for(i in seq(along=dataSource@pathways)){
    p <- dataSource@pathways[i]
    dists[[p]] <- resampleGeneLists(pathway=p,
                                    len=length(geneList), n=n,
                                    dataSource=dataSource)
    setDists[[p]] <- sim2Pathway(pathway=p, doms=listDoms,
                                 dataSource=dataSource)
    ## report progress
    if(verbose && capabilities("tcltk"))
      updateProgress((i)/np*100, autoKill=TRUE,
                     sub= paste("(", i, " of ", np, ")\n",p, sep=""))
  }
  ## compute pvalues from sample distributions
  pvals <- mapply(function(x,y) sum(x>y)/n, dists, setDists)
  pvals[pvals==0] <- 1/n
  res <- list(similarity=setDists, pvalue=pvals)
  if(samples)
      res$dist <- dists
  return(res)
}

