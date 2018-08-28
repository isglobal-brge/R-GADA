\name{getCounts}
\alias{getCounts}

\title{Number of CNVs per sample}
\description{
This function summarize the number of samples that are having an alteration in a given genomic region 
}

\usage{
getCounts(x, group, id="sample")
}

\arguments{
  \item{x}{an object of class 'GRanges' obtained after applying 'getCNVs' function}
  \item{group}{the name of the grouping variable available in the object 'x'}
  \item{id}{the name of the variable having sample id. The default is 'sample'}
}

\value{
This function returns a list having two elements. The first element, called 'counts', contains the number of samples of each segment by the variable 'group'. The second element, called 'n', contains the total number of CNVs by 'group' variable. 
}

\references{
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples". BMC Bioinformatics, 2010;11:380.
}


\seealso{\code{\link{parSBL}}, \code{\link{setupParGADA}}, \code{\link{parBE}}}


\examples{
 \dontrun{
  See the vignette
 } 
}

\keyword{utilities}


