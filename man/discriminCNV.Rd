\name{discriminCNV}
\Rdversion{1.1}
\alias{discriminCNV}
\alias{print.cnv.discrimin}
\alias{plot.cnv.discrimin}

\title{
Discrimination Analysis for CNV data matrices
}
\description{
An ade4 wrapping function to performs a Discriminat Analysis based on multivariate dimensionality reduction such as 
multiple correspondenc analysis, correspondence analysis or principal component analysis. 
Computes the correlation of each variable to the axes defined by each population 
in the principal component subspace. 
}


\usage{

discriminCNV(mat.f, pop.cla, select = NULL, type = "acm", verbose = TRUE, ...)

\method{print}{cnv.discrimin}(x, ...)
\method{plot}{cnv.discrimin}(x, all=FALSE, xax = 1, yax = 2, ...)

}

\arguments{
  \item{mat.f}{\code{data.frame} typically obtained from the matrix reduction performed by \code{\link{reduceMatrix}} or \code{\link{getReducedData}}
}
  \item{pop.cla}{
\code{factors} vector with population labeling of each subject. The subject entries 
must match the entries for the rows of \code{mat.f}  
}
  \item{select}{
list of CNVs to be considered in the discrimination
}
  \item{type}{ for dimensionality reduction
 \code{"acm"} for multiple correspondence analysis
 \code{"coa"} for correspondence analysis
 \code{"pca"} for principal component analysis
}
  \item{verbose}{
whether to print detailed progress on R prompt.
}
  
  \item{x}{ an object of class 'cnv.discrim'}
  \item{all}{similar to 'plot.discrim' of ade4 package. The default is FALSE}
  \item{xax}{number of component plotted in the x-axis. Default is 1st}
  \item{yax}{number of component plotted in the y-axis. Default is 2nd}
  \item{\dots}{
additional parameters for ade4 fucntions
}
}

\value{
Return a list of class \code{cnv.discrimin} with the result of the discrimination 
and the correlation between the CNVs and the directions of all centroid axes.

 \item{proj }{correlation of each CNV onto the centroid directions}
 \item{cen }{centroid directions}
 \item{mat.f }{matrix input}
 \item{pop }{population labeling -input}

}

\references{
Caceres A, Basagana X, Gonzalez JR: Multiple Correspondence Discriminant Analysis: An Application to Detect Stratication in Copy Number Variation 2009. [ISCB Prague, Abstract: S05.03,
submitted to the special issue of Statistics in Medicine].
}

\examples{

\dontrun{

# 'HapMap270reducedData.RData' file can be downloaded from ...

load("HapMap270reducedData.RData")
cp<-discrimin.cnv(mat.f,pop.cla)
plot(cp)

}
}

\keyword{ utilities }

