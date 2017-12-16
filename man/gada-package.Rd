\name{gada-package}
\alias{gada-package}
\alias{gada}
\docType{package}

\title{
Genome Alteration Detection Algorithm (GADA)
}

\description{
This package implements 
a flexible and efficient pipeline to detect copy number alterations from microarray data
(Pique-Regi R and Gonzalez JR, 2008).
The package can import the raw copy number normalized intensities provided by
Illumina BeadStudio, the Affymetrix powertools, or any similar format. 
Probes belonging to different samples are split into separate files and can be analyzed in parallel 
in a cluster/multicore computer. 
The speed and accuracy of the genome alteration detection analysis (GADA) approach 
combined with parallel computing results in one of the fastest and most accurate methods, 
and it is  specially suitable to extract CNA on new genomewide scans involving
hundreds of samples and high density arrays with millions of markers.
The package also provides some functions for summarizing and visualizing 
regions altered at chromosome level (including cytobands).
}

\details{
\tabular{ll}{
Package: \tab gada\cr
Type: \tab Package\cr
Version: \tab 0.7-3\cr
Date: \tab 2008-11-26\cr
License: \tab GPL \cr
URL: \tab http://groups.google.com/group/gadaproject http://www.creal.cat/jrgonzalez/software  \cr
LazyLoad: \tab yes\cr
}

~~ An overview of how to use the package, including the most important ~~
~~ functions ~~

}

\author{
Juan R Gonzalez <jrgonzalez@creal.cat>, Alejandro Caceres <acaceres@creal.cat> and Roger Pique-Regi <rpique@gmail.com>

Maintainer: Juan R Gonzalez <jrgonzalez@creal.cat>

}


\references{
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009

Pique-Regi R, Monso-Varona J,Ortega A, Seeger RC, Triche TJ, Asgharzadeh S. "Sparse representation and Bayesian detection of the genome copy number alterations from microarray data", Bioinformatics , Feb 2008
}

\keyword{ package }

\examples{
\dontrun{
# import data
download.file("http://www.creal.cat/jrgonzalez/GADA/dataIllumina.txt","dataIllumina.txt")

# creating object of class setupGADA
dataIllumina<-setupGADAIllumina(file="dataIllumina.txt", log2ratioCol=5, NumCols=6)

# Segmentation procedure
step1<-SBL(dataIllumina, estim.sigma2=TRUE)
step2<-BackwardElimination(step1,T=4.5,MinSegLen=3) 

# print
step1
step2

# summary
summary(step2)

# plot log2ratio and segments
plotRatio(step2)

# the same for chromosome 12
plotRatio(step2, 12)
}
}
