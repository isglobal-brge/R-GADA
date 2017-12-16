\name{plotWG}
\alias{plotWG}

\title{Plot Whole Genome}

\description{
Plots altered probes, in a given proportion of individuals, for the entire genome when several samples are analyzed
}

\usage{
plotWG(x, min.percentage=0.05, max.number.cnv=100, length.base, chr=c(1:22,"X","Y"))
}

\arguments{
  \item{x}{An object of class 'summaryparGADA'}
  \item{min.percentage}{The minimum percentage of individuals with a shared altered segment to be plotted. The default is 0.05 (e.g., 5\%) }
  \item{max.number.cnv}{For a given chromosome, the maximum number of altered segments found in a given individual. The default is 100. If the number of segments for a given individual is larger than 'max.number.cnv' this sample is discarded.}
  \item{length.base}{a vector with two components. The first corresponds to the minimum size required to consider a segment as a CNV, while the second one corresponds to the maximum size. If missing, this information is obtained from the attribute 'length.base' of object 'x'.}
  \item{chr}{A vector indicating the chromosomes to be plotted. The default is c(1:22, "X", "Y")}
}


\value{
 No return value, just the plot
}

\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}


\seealso{\code{\link{parBE}}, \code{\link{parSBL}}}

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

###
### Third step
### 

# summaryze all samples
allSamples<-summary(myExample)

# plot entire genome
plotWG(allSamples)

# plot a given chromosome
plot(allSamples, chr=6)

}

}

\keyword{utilities}

