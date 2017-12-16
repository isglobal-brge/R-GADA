\name{setupParGADAaffy}
\alias{setupParGADAaffy}

\title{Parallel version of setupGADAaffy}

\description{
This functions calls repeteadly to setupGADAaffy to import log2ratio intensities for several samples to gada
}

\usage{
setupParGADAaffy(folder, files, verbose = TRUE, ...)
}

\arguments{
  \item{folder}{The folder where data is stored. Not required}
  \item{files}{The names of the files with Illumina data. Not required}
  \item{verbose}{Should information about process be printed in the console? The default is TRUE}
  \item{\dots}{ Other arguments passed through '\code{\link{setupGADAaffy}}'}
}

\details{
See setupGADAaffy
}


\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}


\seealso{\code{\link{setupParGADAIllumina}}, \code{\link{parSBL}}, \code{\link{parBE}}, \code{\link{setupGADAaffy}} }

\examples{
\dontrun{
###
### First step (required)
### 

# Read data 
# a folder called rawData containing a file for each individual is required

# Affymetrix data
myExample<-setupParGADAaffy(log2ratioCol=4, NumCols=5)

}

}

\keyword{utilities}

