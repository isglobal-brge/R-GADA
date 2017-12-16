% --- gada-internal.Rd ---

\name{gada-internal}
\alias{gada-internal}
\alias{drawChromosome}
\alias{drawChromosomeZoom}
\alias{plotlogRatio}
\alias{selectSegment}
\alias{collapseInfo}
\alias{WextIextToSegments}
\alias{SBL.fit}
\alias{BackwardElimination.fit}
\alias{plapply}
\alias{countAltered.i}
\alias{findNormalLimits}
\alias{ideogram}
\alias{tamanysChroms}
\alias{get.index}
\alias{get.matrix}
\alias{getProj}
\alias{project.pop}
\alias{assocCNV.i}

\title{Internal gada functions}
\description{Internal gada functions}
\usage{
drawChromosome(chr, size = 0.05, limits, print.names=TRUE, ...)
drawChromosomeZoom(chr, size = 0.05, chromStart, chromEnd, print.names=TRUE, ...)
plotlogRatio(x, chr, chromStart, chromEnd, B.allele.freq=FALSE, num.points, postscript=FALSE, ylim, ...)
selectSegment(x, o)
collapseInfo(x)
WextIextToSegments(x)
SBL.fit(i,x,control,gen.info,chr,sigma2,aAlpha,maxit,tol,debug)
plapply(X, FUN, ...)
countAltered.i(i, x, max.number.cnv, length.base, gen.info,chr.lab)
findNormalLimits(x)
getProj(mat, pop, select=NULL, saveInfo=TRUE, type="acm",verbose=TRUE,...)
project.pop(dd,factor=FALSE,norm=FALSE, verbose=TRUE)
get.matrix(crom.data, probe.pos, verbose,size.min=500,size.max=1.e6)
get.index(pos,probe.pos)
assocCNV.i(x, formula, family=binomial, ...)
}

\details{These are not to be called by the user}
\keyword{internal}


