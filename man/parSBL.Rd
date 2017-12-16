\name{parSBL}
\alias{parSBL}

\title{Parallel SBL segmentation}

\description{
This function perform SBL segmentation for a collection of samples. 
}

\usage{
parSBL(x, Samples, estim.sigma2, aAlpha, verbose = TRUE, ...)
}

\arguments{
  \item{x}{an object of class 'parGADA' obtained using 'setupParGADAIllumina' or 'setupParGADAaffy' functions}
  \item{Samples}{the samples to be analyzed}
  \item{estim.sigma2}{array noise level (see 'SBL' function)}
  \item{aAlpha}{sparseness hyperparameter (see 'SBL' function)}
  \item{verbose}{is information about each individual showed in the console? Default is TRUE}
  \item{\dots}{other arguments passed through the function 'SBL'}
}

\details{
This function calls repeteadly to function 'SBL' to obtain all segments for each sample. It can be parallelized when multiple processors are available. This has been implemented using 'snow' package: http://www.sfu.ca/sblay/R/snow.html. See vignettes for further details 
}

\value{
This function saves the segments for each sample in a directory called SBL
}

\references{
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}

\seealso{\code{\link{parBE}}, \code{\link{setupParGADAIllumina}}, \code{\link{setupParGADAaffy}}}


\examples{
\dontrun{
###
### First step (required)
### 

# Read data 
# a folder called rawData containing a file for each individual is required

# Illumina data
myExample<-setupParGADAIllumina(log2ratioCol=4, NumCols=5)

# Affymetrix data
myExample<-setupParGADAaffy(log2ratioCol=4, NumCols=5)

###
### Second step
### 

# segmentation for all samples
parSBL(myExample, estim.sigma2=TRUE, aAlpha=0.8)
parBE(myExample,T=8, MinSegLen=8)
}

}

\keyword{utilities}

