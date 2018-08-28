\name{setupParGADA}
\alias{setupParGADA}

\title{Parallel version of setupGADA}

\description{
This functions calls repeteadly to setupGADA to import raw data (pennCNV format) for several samples to R-GADA
}

\usage{
setupParGADA(folder, files, MarkerIdCol=1,
             ChrNameCol=2, ChrPosCol=3,
             chrs=c(1:22, "X", "Y"),
             XY=TRUE, verbose=TRUE, 
             mc.cores=1, ...)

}

\arguments{
  \item{folder}{The folder where data is stored. Not required}
  \item{files}{The names of the files with Illumina data. Not required}
  \item{MarkerIdCol}{The column having marker ids in the pennCNV files. The defaul is 1}
  \item{ChrNameCol}{The column having chromosome names in the pennCNV files. The defaul is 2}
  \item{ChrPosCol}{The column having genomic position in the pennCNV files. The defaul is 3}
  \item{chrs}{A vector containg the available chromosomes. The default is c(1:22, "X", "Y")}
  \item{XY}{Do data have X and Y labels in annotation (instead of 23 and 24)? The default is TRUE}
  \item{verbose}{Should information about process be printed in the console? The default is TRUE}
  \item{mc.cores}{number of cores to be used when using multiple cores (see argument 'mc.cores' from 'mclapply' function of 'parallel' library) }
  \item{\dots}{ Other arguments passed through 'setupGADA'}
}

\details{
See setupGADA
}


\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples". BMC Bioinformatics, 2010;11:380.
}

\seealso{\code{\link{setupParGADA}}, \code{\link{parSBL}}, \code{\link{parBE}} }

\examples{
 \dontrun{
  See the vignette
 } 
}

\keyword{utilities}

