\name{BackwardElimination}
\alias{BackwardElimination}
\alias{print.BackwardElimination}
\alias{summary.BackwardElimination}

\title{Backward Elimination of segmentation breakpoints}

\description{
  Adjusts the breakpoint sensitivity vs. FDR trade-off after fitting the
  SBL model
}


\usage{
BackwardElimination(x, T, MinSegLen, saveInfo = TRUE)
}


\arguments{
  \item{x}{an object of class 'SBL'}
  \item{T}{critical value of the backward elimination algorithm}
  \item{MinSegLen}{minimum number of probes each segment must contain}
  \item{saveInfo}{Indicates whether annotation data is kept. The default is TRUE}
}

\details{

  After fitting the SBL model, we can use this function to quickly
  adjust the FDR.

  The 'T' argument indicates the critical value of the BE algorithm. That is,
  the statistical score $t_m$ associated with each breakpoint $m$ remaining in the model
  has to be higher than $T$. The score $t_m$ can be interpreted as the
  difference between the sample averages of the probes falling on the left and right
  segment, divided by a pooled estimation of the standard error.
  Asympotically, when the number of probes on the right and left segments are
  very large this score will converge to a standard normal distribution,
  i.e. N(0,1).
  
  The argument 'MinSegLen' can be used to limit the minimum number
  of probes each CNA segment must contain. We recommend setting at least
  'MinSegLen=3' (default value) to eliminate false detections due to
  extreme outliers. 

}

\value{
   An object of class 'BackwardElimination'

  'print' returns the number of discontinuities or segments by chromosome 
  
  'summary' returns the altered segments (gains and loses) with the begining and the end of each one, as well as the number of probes in each segment. This information is given by chromosome

  'plotRatio' produces a plot of log2ratios and segments by chromosome. See \code{\link{plotRatio}} for
  further details.

  See examples for further illustration about all previous issues. 


}

\references{ 

Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009

Pique-Regi R, Monso-Varona J,Ortega A, Seeger RC, Triche TJ, Asgharzadeh S. "Sparse representation and Bayesian detection of the genome copy number alterations from microarray data", Bioinformatics , Feb 2008

}


\seealso{\code{\link{setupGADAIllumina}}, \code{\link{setupGADAaffy}}, \code{\link{setupGADAgeneral}}, \code{\link{SBL}}, \code{\link{parSBL}}, \code{\link{parBE}}}


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

\keyword{models}
