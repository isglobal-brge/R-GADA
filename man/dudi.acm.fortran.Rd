\name{dudi.acm.fortran}
\Rdversion{1.1}
\alias{dudi.acm.fortran}
\title{
Duality diagram for Multiple Correspondenc Analysis
}
\description{
optimized \code{dudi.acm} from the ade4 package. The subfunction \code{floc} has
been reprogramed in FORTRAN, to increase speed.
}
\usage{
dudi.acm.fortran(df, row.w = rep(1, nrow(df)), scannf = TRUE, nf = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
\code{data.frame} containing only factors
}
  \item{row.w}{
vector of row weights, by default uniform weighting
}
  \item{scannf}{
a logical value indicating whether the eigenvalues bar plot should be displayed
}
  \item{nf}{
if sacnnf FALSE, and integer indicating number of kept axes
}
}
\details{
see ade4 documentation for further information
}

\value{
returs a list of class \code{acm} and \code{dudi}
}


\references{
see ade4 website: http://pbil.univ-lyon1.fr/ADE-4/
}

\examples{

\dontrun{

# 'HapMap270reducedData.RData' file can be downloaded from ...

load("HapMap270reducedData.RData")
cm<-dudi.acm.fortran(mat.f,scan=FALSE)
plot(cm, pop.cla, var=FALSE,pnt=0.7)

}

}

\keyword{ utilities }
