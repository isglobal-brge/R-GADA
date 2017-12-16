\name{matrixPlot}
\Rdversion{1.1}
\alias{matrixPlot}

\title{
Matrix image of CNV blocks
}

\description{
Plots the matrix of CNV blocks }

\usage{
matrixPlot(x, fac = TRUE, ...)
}

\arguments{
  \item{x}{
\code{matrix} or \code{data.frame} of CNV variables, coding three levels each:
-1:loss, 0:no-change, 1:gain in copy number variation.
}
  \item{fac}{
whether the x columns encode factors or numeric values
}
  \item{\dots}{
additional parameters than can be passed to \code{matrixPlot} 
}
}

\examples{

\dontrun{

# 'Hapmap270Cnv.RData' file can be downloaded from ...

load("Hapmap270Cnv.RData")

ch<-6
mat<-get.chromosome.dat(segments,ch,gen.info)
mat.f<-reduce.matrix(mat,gen.info,ch,varSimil=0.99,subVariation=0.99)
cnv.blocks<-attr(mat.f,"cnv.blocks")
matrixPlot(mat.f)
}

}

\keyword{ utilities }

