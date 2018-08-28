\name{plotCNVs}
\alias{plotCNVs}

\title{Gviz plot of CNVs}
\description{
This function plots the CNVs of each sample in a given genomic region and annotates the genes using Gviz 
}

\usage{
plotCNVs(x, range, genome="hg19", drawGenes=FALSE, group,  
         mosaic = FALSE)
}

\arguments{
  \item{x}{an object of class 'cnvGRanges' obtained after applying 'getCNVs' function}
  \item{range}{a 'GRanges' object indicating the genomic region to be depicted}
  \item{genome}{the reference genome. The default is 'hg19'}
  \item{drawGenes}{should genes be showed in the plot? The default is FALSE}
  \item{group}{a grouping variable to provide CNVs for each level. Only to levels are allowed. Not required.}
  \item{mosaic}{argument not relevant for R-GADA (see MAD package)}
}

\value{
This function returns an Gviz plot showing the CNVs of each sample in a given genomic region. Gains are in green color, while loses are printed in blue.
}

\references{
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples". BMC Bioinformatics, 2010;11:380.
}


\seealso{\code{\link{getCNVs}}}


\examples{
 \dontrun{
  See the vignette
 } 
}

\keyword{utilities}


