\name{setupParGADAIllumina}
\alias{setupParGADAIllumina}

\title{Parallel version of setupGADAIllumina}

\description{
This functions calls repeteadly to setupGADAIllumina to import log2ratio intensities for several samples to gada
}

\usage{
setupParGADAIllumina(folder, files, XY=TRUE, verbose = TRUE, sort=TRUE, ...)
}

\arguments{
  \item{folder}{The folder where data is stored. Not required}
  \item{files}{The names of the files with Illumina data. Not required}
  \item{XY}{Do data have X and Y labels in annotation or 23 and 24? The default is TRUE}
  \item{verbose}{Should information about process be printed in the console? The default is TRUE}
  \item{sort}{Should data be sorted by genomic position? See argument 'sort' in setupGADAIllumina'}
  \item{\dots}{ Other arguments passed through 'setupGADAIllumina'}
}

\details{
See setupGADAIllumina
}


\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}

\seealso{\code{\link{setupParGADAIllumina}}, \code{\link{parSBL}}, \code{\link{parBE}} }

\examples{
\dontrun{
###
### First step (required)
### 

# Read data 
# a folder called rawData containing a file for each individual is required

# Illumina data
myExample<-setupParGADAIllumina(log2ratioCol=4, NumCols=5)
}

}

\keyword{utilities}

