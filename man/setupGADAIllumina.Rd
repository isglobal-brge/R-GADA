\name{setupGADAIllumina}
\alias{setupGADAIllumina}

\title{Imports Illumina array data to gada}

\description{
This function encapsulates information obtained from Illumina platform (BeadStudio format) to be analyzed with gada
}

\usage{
setupGADAIllumina(file, NumCols, log2ratioCol, MarkerIdCol = 1, ChrNameCol = 2, ChrPosCol = 3, XY=TRUE, sort = TRUE, orderProbes, saveGenInfo = TRUE)
}

\arguments{
  \item{file}{A file with Illumina data (obtained from BeadStudio) }
  \item{NumCols}{The number of columns in 'file'}
  \item{log2ratioCol}{The column in 'file' containing log2ratio intensities}
  \item{MarkerIdCol}{The column in 'file' containing the name of the marker}
  \item{ChrNameCol}{The column in 'file' containing the chromosome}
  \item{ChrPosCol}{The column in 'file' containing the genomic position}
  \item{XY}{Do data have X and Y labels in annotation or 23 and 24? The default is TRUE}
  \item{sort}{Should data be sorted using annotation data? The default is TRUE}
  \item{orderProbes}{Used when sort is TRUE. See details}
  \item{saveGenInfo}{Should annotation data be saved? The default is TRUE. See details}
}

\details{
The arguments 'orderProbes' and 'saveGenInfo' should not be changed by the user.
}

\value{
An object of class 'setupGADA'
}

\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}

\seealso{ \code{\link{setupGADAgeneral}}, \code{\link{setupGADAaffy}}}

\examples{
\dontrun{
# get data
download.file("http://www.creal.cat/jrgonzalez/GADA/dataIllumina.txt",
               "./dataIllumina.txt")

dataIllumina<-setupGADAIllumina(file="dataIllumina.txt",log2ratioCol=5,NumCols=6)
dataIllumina
}
}

\keyword{utilities}

