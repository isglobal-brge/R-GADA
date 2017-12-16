\name{exportToBED}
\alias{exportToBED}

\title{Export segments to BED format}

\description{
Export information about altered segments to BED format

}

\usage{
exportToBED(x, num=300, strand="+", col.gain="0,255,255", col.loss="255,0,0", file="BED.txt")
}

\arguments{
  \item{x}{An object of class 'summaryparGADA'}
  \item{num}{}
  \item{strand}{}
  \item{col.gain}{}
  \item{col.loss}{}
  \item{file}{The name of the file to be saved. The default is 'BED.txt'}
}


\value{
 No return value, just a file 
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

# export to BED
exportToBED(allSamples)

}

}

\keyword{utilities}

