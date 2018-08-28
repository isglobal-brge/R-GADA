\name{getCNVs}
\alias{getCNVs}

\title{CNVs after performing calling with R-GADA}
\description{
This function provides the CNVs of each sample by chromosome after performing segmentation procedure using GADA. 
}

\usage{
getCNVs(x)
}

\arguments{
  \item{x}{an object of class 'parGADA' obtained after applying 'parBE' function}
}

\value{
This function returns an object of class 'cnvGR' which is a 'GRanges' object that contains all the segments detected as CNV gain or loss.
}

\references{
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples". BMC Bioinformatics, 2010;11:380.
}


\seealso{\code{\link{parSBL}}, \code{\link{setupParGADA}}, \code{\link{parBE}}
}

\examples{
 \dontrun{
  See the vignette
 } 
}

\keyword{utilities}


