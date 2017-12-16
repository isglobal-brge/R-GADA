\name{getAlteredProbes}
\alias{getAlteredProbes}

\title{Get Altered Probes}

\description{
This function returns a list with two components. Each component contains a data frame with the name of the probe and the number of individuals who have this probe altered. First component contains such information for gains and second one for losses. 

}

\usage{
getAlteredProbes(x, chr, min.perc = 0.1, max.number.cnv = 100, length.base)
}

\arguments{
  \item{x}{an object of class 'summaryparGADA'}
  \item{chr}{the desired chromosome}
  \item{min.perc}{minimum percentage of individuals that must have a given altered probe. By default is 0.1 (e.g. 10\%)}
  \item{max.number.cnv}{maximum number of segments for each individual. The default is 100. If the number of segments for a given individual is larger than 'max.number.cnv' this sample is discarded}
  \item{length.base}{a vector with two components. The first corresponds to the minimum size required to consider a segment as a CNV, while the second one corresponds to the maximum size. If missing, this information is obtained from the attribute 'length.base' of object 'x'}
}


\value{
  A list with the following components

  \item{gains}{a data frame containning the name of probes with a gain and the number of individuals who have this probe altered}
  \item{losses}{a data frame containning the name of probes with a loss and the number of individuals who have this probe altered}
}

\references{
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}


\seealso{ \code{\link{parBE}}, \code{\link{parSBL}}}

\examples{
# see vignette
}

\keyword{utilities}
