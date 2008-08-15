\name{gseDomain}
\alias{gseDomain}

\title{Geneset enrichment based on InterPro domain signatures}
\description{
  Compute the similarity to pathways specified trough \code{dataSource}
  for a set of entrezgene identifiers.
}
\usage{
gseDomain(dataSource, geneset, n=10000, verbose=TRUE, samples=FALSE)
}
\arguments{
  \item{dataSource}{Object of class \code{ipDataSource} containing
pathway and InterPro domain mappings}
  \item{geneset}{Character vector of entrezgene identifiers}
  \item{n}{Number of subsampling iterations}
  \item{verbose}{Toggle progress report}
  \item{samples}{Logical indicating whether to return the similarity
    measures for all the resamples.}
}
\details{
 Use this function to compute p-values for similarity of the domain
signature of a gene set to all signatures of the pathways defined in
\code{dataSource}. You should have created \code{dataSource} using
either function \code{\link{dataSource}} or \code{\link{getKEGGdata}}.
}
\value{
  A list with items
  
  \item{similarity }{Named vector of similarity measures for each
    pathway}
 
  \item{pvalue }{The p-values of similarity to each pathway. A named
    vector.}
 
  and optional item (if \code{samples=TRUE})
  
  \item{dist }{A named list containing similarity measures for all the
    resamples}
}

\author{Florian Hahne}

\seealso{\code{\link{gseDomain}}}
\examples{
## see Vignette of this package for examples how to use this function
}
\keyword{manip}