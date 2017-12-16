\name{rankVariables}
\Rdversion{1.1}
\alias{rankVariables}
\title{
CNV ranking
}

\description{
Ranks each CNV block according to their highest correlation across population centroids. 
Hence, it assigns the population that the variable discriminate the best, according to 
the discriminant analysis performed by \code{\link{discriminCNV}}.
}

\usage{
rankVariables(x, pop, cnv.blocks = NULL, proj = TRUE, ...)
}

\arguments{
  \item{x}{
object of class \code{cnv.discrimin} or reduced matrix as obtained by \code{\link{reduceMatrix}} or \code{\link{getReducedData}}.
}
  \item{pop}{
\code{factors} vector with population labeling of each subject. 
}
  \item{cnv.blocks}{
as obtained by the attribute "cnv.blocks" from the result of
\code{\link{reduceMatrix}} or \code{\link{getReducedData}}.
}
  \item{proj}{
whether projection on centroid direction or principal axes should be used for the ranking.
}
  \item{\dots}{
additional parameters passed to \code{\link{discriminCNV}}
}
}

\value{
\code{rankVariables} produces a \code{data.frame} in which each CNV block is 
listed with the number of CNV is composed of, the initial and final genomic possition, 
the chromosome it belongs to, and the correlation to the axes defined by its most 
prominent population.   
}
\references{
Caceres A, Basaga~na X, Gonzalez JR: Multiple Correspondence Discriminant Analysis: An
Application to Detect Stratication in Copy Number Variation 2009. [ISCB Prague, Abstract: S05.03,
submitted to the special issue of Statistics in Medicine].
}
\author{
Alejandro Cacereres
}
\seealso{
\code{\link{discriminCNV}}
}
\examples{


\dontrun{

# 'HapMap270reducedData.RData' file can be downloaded from ...

load("HapMap270reducedData.RData")

cp<-discrimin.cnv(mat.f,pop.cla)

var.proj<-rank.variables(cp,cnv.blocks=cnv.blocks)
var.proj[1:10,]
}

}
\keyword{ utilities }

