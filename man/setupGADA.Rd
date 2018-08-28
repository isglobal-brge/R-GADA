\name{setupGADA}
\alias{setupGADA}

\title{Imports array data (pennCNV format) to R-GADA}

\description{
This function encapsulates information in pennCNV format to be analyzed with R-GADA
}

\usage{
setupGADA(file, MarkerIdCol=1, ChrNameCol=2, ChrPosCol=3,
          log2ratioCol, BAFcol, chrs = c(1:22, "X", "Y"),
          XY=TRUE, sort = TRUE, orderProbes, gen.info, 
          saveGenInfo = TRUE) 
}

\arguments{
  \item{file}{A file with Illumina data (obtained from BeadStudio) }
  \item{MarkerIdCol}{The column in 'file' containing the name of the marker}
  \item{ChrNameCol}{The column in 'file' containing the chromosome}
  \item{ChrPosCol}{The column in 'file' containing the genomic position}
  \item{log2ratioCol}{The column in 'file' containing log2ratio intensities}
  \item{BAFcol}{The column in 'file' containing B-allele freqencies}
  \item{chrs}{A vector containg the available chromosomes. The default is c(1:22, "X", "Y")}
  \item{XY}{Do data have X and Y labels in annotation (instead of 23 and 24)? The default is TRUE}
  \item{sort}{Should data be sorted using genomic information? The default is TRUE}
  \item{orderProbes}{Used to avoid sorting}
  \item{gen.info}{genomic information. Not required}
  \item{saveGenInfo}{Should annotation data be saved? The default is TRUE. See details}
}

\details{
The arguments 'orderProbes' and 'saveGenInfo' should not be changed by the user.
}

\value{
An object of class 'setupGADA'
}

\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples". BMC Bioinformatics, 2010;11:380.
}

\examples{
 \dontrun{
  See the vignette
 } 
}

\keyword{utilities}

