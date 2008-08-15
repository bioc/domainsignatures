## Compute the chisquare-based similarity of a gene list to
## a pathway

sim2Pathway <- function(pathway, doms, dataSource){
  ## the unique pathway domains
  ps <- dataSource@path2Domains[[pathway]]
  lps <- length(ps)
  ## number of unique domains in universe (cumsum of cont table)
  lu <- dataSource@dims[3]
  ## the contingency table
  x <- sum(doms %in% ps)
  y <- lps-x
  n1 <- length(doms)
  n2 <- lu - n1
  ## modified chisquare test  
  pval <- sage.test(x,y,n1,n2)
  return(-log10(pval))
}

