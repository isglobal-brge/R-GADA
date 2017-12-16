\name{splitDataBeadStudio}
\alias{splitDataBeadStudio}

\title{Split BeadStudio data}

\description{
This function split a file obtained from BeadStudio into different files corresponding to each individual
}

\usage{
splitDataBeadStudio(file, Samples, NumCols, verbose = TRUE)
}

\arguments{
  \item{file}{An unique BeadStudio file containing information for all samples}
  \item{Samples}{The number of individuals}
  \item{NumCols}{The number of Columns for each individual}
  \item{verbose}{Should information about process be printed in the console? The default is TRUE}
}


\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}

\seealso{\code{\link{setupGADAIllumina}}}

\examples{
\dontrun{
# get data
download.file("http://www.creal.cat/jrgonzalez/GADA/exampleBeadStudio.txt","./exampleBeadStudio.txt")

splitDataBeadStudio("exampleBeadStudio.txt",Samples=8,NumCols=5)

}

}


\keyword{utilities}

