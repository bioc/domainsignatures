\name{getKEGGdescription}
\alias{getKEGGdescription}

\title{KEGG description from ID}
\description{
 Get description of KEGG pathways from a list of KEGG identifiers.
}
\usage{
getKEGGdescription(ids)
}

\arguments{
  \item{ids}{Character vector of KEGG identifiers}
}

\value{
  Character vector of KEGG descriptions
}

\author{Florian Hahne}

\examples{
getKEGGdescription("hsa03050")
}
\keyword{IO}