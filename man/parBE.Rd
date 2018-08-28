\name{parBE}
\alias{parBE}
\alias{summary.parGADA}
\alias{labels.parGADA}
\alias{labels.summaryParGADA}
\alias{print.summaryParGADA}

\title{Parallel Backward Elimination}
\description{
This function performs backward elimination procedure after SBL segmentation for a collection of samples. 
}

\usage{
parBE(x, Samples, T, MinSegLen, verbose = TRUE, 
      mc.cores=1, ...)

\method{print}{summaryParGADA}(x, ...)
\method{labels}{parGADA}(object, ...)
\method{labels}{summaryParGADA}(object, ...)

}

\arguments{
  \item{x}{an object of class 'parGADA' obtained using 'setupParGADAIllumina' or 'setupParGADAaffy' functions}
  \item{Samples}{the samples to be analyzed}
  \item{T}{critical value of the backward elimination algorithm (see 'BackwardElimination' function)}
  \item{MinSegLen}{minimum number of probes each segment must contain (see 'BackwardElimination' function)}
  \item{verbose}{is information about each individual showed in the console? Default is TRUE}
  \item{object}{an object of class 'parGADA' or 'summaryParGADA'}
  \item{mc.cores}{number of cores to be used when using multiple cores (see argument 'mc.cores' from 'mclapply' function of 'parallel' library) }
  \item{\dots}{other arguments passed through the function 'BackwardElimination'}
}

\details{
This function calls repeteadly to function 'BackwardElimination' after being used the function 'parSBL' to obtain all segments for each sample. It can be parallelized when multiple processors are available. This has been implemented using 'snow' package: http://www.sfu.ca/sblay/R/snow.html. See vignettes for further details 
}

\value{
This function saves the segments for each sample in a directory called SBL
}

\references{
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples". BMC Bioinformatics, 2010;11:380.
}


\seealso{\code{\link{parSBL}}, \code{\link{setupParGADA}}}

\examples{
 \dontrun{
  See the vignette
 } 
}

\keyword{utilities}


