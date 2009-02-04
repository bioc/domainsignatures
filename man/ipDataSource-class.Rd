\name{ipDataSource-class}
\docType{class}
\alias{ipDataSource-class}
\alias{ipDataSource}
\alias{show,ipDataSource-method}



\title{A class to store mapping information between genes, pathways and
interPro domains}

\description{This class represents the data necessary to run
\code{\link{gseDomain}}}


\section{Creating Objects}{
Objects can be created using\cr
\code{  new('ipDataSource',}\cr
\code{    genes        = ...., # Object of class character}\cr
\code{    pathways     = ....  # Object of class character}\cr
\code{    domains      = ....  # Object of class character}\cr
\code{    gene2Domains = ....  # Object of class environment}\cr
\code{    path2Domains = ....  # Object of class environment}\cr
\code{    type         = ...., # Object of class character}\cr
\code{  )}\cr
or the functions \code{\link{dataSource}} or \code{\link{getKEGGdata}}.}

\section{Slots}{\describe{
\item{\code{genes}:}{Vector of unique entrezgene identifiers}
\item{\code{pathways}:}{Vector of unique pathway identifiers}
\item{\code{domains}:}{Vector of unique InterPro identifiers}
\item{\code{gene2Domains}:}{Hash table mapping entrezgene IDs to
  Interpro IDs}
\item{\code{path2Domains}:}{Hash table mapping pathway IDs to InterPro
  IDs}
\item{\code{type}:}{Type of pathway. A character skalar}
}}

\section{Methods}{\describe{
\item{show}{display summary.}
}}

\details{
You should always create these objects using either one of the functions \code{\link{dataSource}} or \code{\link{getKEGGdata}}.
}
			
\author{Florian Hahne}

\seealso{\code{\link{dataSource}} or \code{\link{getKEGGdata}}
}


\keyword{classes}

