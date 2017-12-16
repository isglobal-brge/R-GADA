\name{getNamesProbes}
\Rdversion{1.1}
\alias{getNamesProbes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CNVs with highes correlations
}

\description{
Obtains the CNV blocks that have the highes correlation according to 
\code{rankVariables}.
}

\usage{
getNamesProbes(x, min.correlation = 0.3)
}

\arguments{
  \item{x}{  \code{data.frame} with the ranking of the CNV blocks according to their correlation  to the population centroids. } 
  \item{min.correlation}{value of the minimum correlation to select the variables from the ranking.}
}


\details{
The output of \code{getNamesProbes} can be used to re-train a discriminant analysis \code{discriminCNV}. 
}


\value{
produces a vector of characters with the CNV blocks that are most relevant in the discrimination summarized on \code{x}.
}


\examples{

\dontrun{

# 'HapMap270reducedData.RData' file can be downloaded from ...

load("HapMap270reducedData.RData")

cp<-discrimin.cnv(mat.f,pop.cla)

var.proj<-rank.variables(cp,cnv.blocks=cnv.blocks)
select<-getNamesProbes(var.proj, min.correlation=0.5)
#re-train discriminant
cp.2<-discrimin.cnv(mat.f, pop.cla, select=select)
plot(cp.2)

#re-compute correspondence analysis with selected variables.
cm<-dudi.acm.fortran(mat.f[select], scan=FALSE, nf=3)
plot(cm)
}
}
\keyword{ utilities }

