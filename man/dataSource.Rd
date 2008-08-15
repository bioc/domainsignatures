\name{dataSource}
\alias{dataSource}

\title{Contructor for annotation object}
\description{
  This function creates the necessary annotation object of class
  \code{ipDataSource} containing the pathway and InterPro ID
  mappings.
}
\usage{
dataSource(mapping, type = "generic")
}

\arguments{
  \item{mapping}{A named list providing a mapping between enterzgene
    identifiers and arbitrary groupings of genes or pathways.}
  \item{type}{The type of pathway. A character skalar.}
}
\details{
  For genes without pathway membership, \code{NA} list items need to be
  included in \code{mapping}. The names of the list comprise the gene
  universe to test against. The function will access the ensembl
  \code{biomaRt} database in order to retrieve the necessary InterPro domain
  information.
}
\value{
  Object of class \code{ipDataSource}
}
\author{Florian Hahne}

\seealso{\code{\link{gseDomain}}, \code{\link{getKEGGdata}}}
\examples{
grouping <- list("653361"=c("pw1", "pw2"), "729230"="pw1",
                 "415117"="pw3") 
dataSource(grouping)
}
\keyword{IO}