\name{plotRatio}
\alias{plotRatio}
\alias{plotRatio.setupGADA}
\alias{plotRatio.BackwardElimination}
\alias{plotRatio.parGADA}


\title{Plot log2Ratio intensities}

\description{
This function produces a plot of log2ratios and segments by chromosome
}

\usage{
plotRatio(x, \dots)

\method{plotRatio}{setupGADA}(x, chr, chromStart, chromEnd, B.allele.freq=FALSE, \dots)

\method{plotRatio}{BackwardElimination}(x, chr, num.points, postscript=FALSE, \dots)

\method{plotRatio}{parGADA}(x, Sample, segments=FALSE, \dots)

}


\arguments{
  \item{x}{An object of class 'setupGADA', 'BackwardElimination', or 'parGADA'}

  \item{chr}{A given cromosome}
  \item{chromStart}{The beginning of the cromosome}
  \item{chromEnd}{The end of the cromosome}
  \item{B.allele.freq}{The function is prepared to show B-allele frequency for Illumina data in future versions. This should not be changed. The defaul is FALSE}
  \item{num.points}{The desired number of log2ratio intensities to be showed in the plot. This is useful for importing figures}
  \item{postscript}{Is the plot exported to a encapsulated postscript file? the default is FALSE}

  \item{Sample}{The sample to be plotted. This is required when an object of class 'parGADA' is passed through the argument 'x'}

  \item{segments}{Are segments added to the plot? This can be used when an object of class 'parGADA' is passed through the argument 'x'}

  \item{\dots}{Other arguments passes through the function plotlogRatio}
}


\value{
 No return value, just the plot
}


\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}

\seealso{\code{\link{setupGADAIllumina}}, \code{\link{setupGADAaffy}}, \code{\link{SBL}}, \code{\link{BackwardElimination}} }

\examples{
\dontrun{
## Simulated data
set.seed(123456)
cn<-rep(c(rep(1,1E5-100),rep(1,100),rep(1,1E5)),4) #Underlying copy number
arrayData<-rnorm(length(cn),mean=(log2(cn)-1),sd=1) #Simulated array
dataSim<-setupGADAgeneral(arrayData)
plotRatio(dataSim)

## Real data
# import data
download.file("http://www.creal.cat/jrgonzalez/GADA/dataIllumina.txt","dataIllumina.txt")

# creating object of class setupGADA
dataIllumina<-setupGADAIllumina(file="dataIllumina.txt", log2ratioCol=5, NumCols=6)

# Segmentation procedure
step1<-SBL(dataIllumina, estim.sigma2=TRUE)
step2<-BackwardElimination(step1,T=4.5,MinSegLen=3) 

# plot log2ratio intensities
plotRatio(dataIllumina)

# plot log2ratio intensities for chromosome 12
plotRatio(dataIllumina, chr=12)

# plot log2ratio intensities and segments
plotRatio(step2)
}

}

\keyword{utilities}

