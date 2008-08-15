\name{getKEGGdata}
\alias{getKEGGdata}

\title{Fetch KEGG annotations and InterPro domains}
\description{
  Get all available KEGG annotations and InterPro domains for a set of
  entrezgene identifiers from the \code{KEGG} annotation package and
  from the ensembl \code{biomaRt}.
}
\usage{
getKEGGdata(universe, pathways)
}

\arguments{
  
  \item{universe}{Character vector of entrezgene identifiers. This is
    the global universe of genes to test against.}

  \item{pathways}{Optional character vector of KEGG pathway
    identifiers. This can be used in order to test for
    over-representation of only a subset of all the available KEGG
    pathways. }

}


\details{ This function is a wrapper around the \code{KEGG} annotation
  package and a customized query of the ensembl \code{biomaRt}
  database. For the gene identifiers in \code{universe} and all KEGG
  pathways it will fetch the necessary information and create an object
  of class \code{ipDataSource} which can later be used as input to
  \code{\link{gseDomain}} }

\value{
  An object of class \code{ipDataSource}.
}

\author{Florian Hahne}

\seealso{\code{\link{gseDomain}}}
\examples{
getKEGGdata()
}
\keyword{IO}