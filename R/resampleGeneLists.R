## Resample n random geneLists of length 'len' from all
## universe genes and compute chisquare-based distance measures
## to signature of 'pathway'

resampleGeneLists <- function(pathway, len, n, dataSource){
  ## the unique pathway domains
  #ps <- get(pathway, dataSource@path2Domains)
  ps <- dataSource@path2Domains[[pathway]]
  lps <- length(ps)
  ## number of unique domains in universe (cumsum of cont table)
  lu <- length(dataSource@domains)
  sGenes <- replicate(n, sample(seq_along(dataSource@genes), len), 
                      simplify=FALSE)
  gs <- lapply(sGenes, function(x,y) unique(unlist(y[x], use.names=FALSE)),
               dataSource@gene2Domains)
  x <- sapply(gs, function(x,y) sum(x%in% y), ps)
  y <- lps-x
  n1 <- listLen(gs)
  n2 <- lu-n1
  pval <- sage.test(x,y,n1,n2)
  return(-log10(pval))
}
