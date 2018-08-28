\name{parSBL}
\alias{parSBL}

\title{Parallel SBL segmentation}

\description{
This function perform SBL segmentation for a collection of samples. 
}

\usage{
parSBL(x, Samples, estim.sigma2, aAlpha, verbose = TRUE, 
       mc.cores=1,...)
}

\arguments{
  \item{x}{an object of class 'parGADA' obtained using 'setupParGADAIllumina' or 'setupParGADAaffy' functions}
  \item{Samples}{the samples to be analyzed}
  \item{estim.sigma2}{array noise level (see 'SBL' function)}
  \item{aAlpha}{sparseness hyperparameter (see 'SBL' function)}
  \item{verbose}{is information about each individual showed in the console? Default is TRUE}
  \item{mc.cores}{number of cores to be used when using multiple cores (see argument 'mc.cores' from 'mclapply' function of 'parallel' library) }
  \item{\dots}{other arguments passed through the function 'SBL'}
}

\details{
This function calls repeteadly to function 'SBL' to obtain all segments for each sample. It can be parallelized when multiple processors are available. This has been implemented using 'snow' package: http://www.sfu.ca/sblay/R/snow.html. See vignettes for further details 
}

\value{
This function saves the segments for each sample in a directory called SBL
}

\references{
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples". BMC Bioinformatics, 2010;11:380.
}

\seealso{\code{\link{parBE}}, \code{\link{setupParGADA}}, \code{\link{setupParGADA}}}


\examples{
 \dontrun{
  See the vignette
 } 
}

\keyword{utilities}

