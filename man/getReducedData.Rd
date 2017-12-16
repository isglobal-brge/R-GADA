\name{getReducedData}
\Rdversion{1.1}
\alias{getReducedData}
\alias{sortPopData}


\title{
Reduced matrix of CNV blocks across the whole genome
}
\description{
Obtains the reduced matrix of CNV segments and subject samples from segmentation 
output. The initail matrix is obtained and reduced for each cromosome in turn. 
The resulting matrix is a concatenation of the reduced matrices for each chromosome.
}
\usage{
getReducedData(x, gg, varSimil = 0.99, subVariation = 0.99)
}
\arguments{
  \item{x}{
concatenated list containing segmentation outcome for each subject, chromosome 
and population.  
}
  \item{gg}{
genetic info: \code{data.frame} with probe names, chromosome number and 
genomic positions of each CNV. 
}
 \item{varSimil}{
\code{varSimil=0.99} produces blocks of neighboring CNVs which do differ 
in more than 99\% across the subject sample. 
}
  \item{subVariation}{
\code{subVariation = 0.99} gets rid of CNV blocks that are constant across the 
the whole subject sample, except for 1\% of the subjects.
}

}
\details{
For \code{x}, concatenation by subject and chromosome, within a single population,
is readily provided by the function \code{summary.parGADA} on an object \code{parGADA} 
that has been previously segemnted. Population concatenation is the result of   
\code{c(summary{parGADA1}), summary{parGADA2}), ...} across groups. Concatenation can be obtained with the function sortPopData(segments,nchr,npop), where segments is the output of GADA, nchr the number of chromosomes analyzed and npop the number of populations in the data set.

}

\value{
\item{getReducedData }{ produces a \code{data.frame} with colums \code{as.factors} representing each CNV block}
\item{attr(getReducedData,"cnv.blocks"}{lists the genetic information of each CNV block}
}


\author{
Alejandro Caceres}
\seealso{
\code{\link{reduceMatrix}}, \code{\link{getChromosomeDat}}
}
\examples{

\dontrun{

# 'Hapmap270Cnv.RData' file can be downloaded from ...

load("Hapmap270Cnv.RData")

mat.f<-getReducedData(segments, gen.info, varSimil=0.99, subVariation=0.90)
cnv.blocks<-attr(mat.f,"cnv.blocks")
matrixPlot(mat.f)

#if segments are available form gada output first sort data by populations then get reduced data.
#segments<-sortPopData(segments,nchr,npop)
#mat.f<-getReducedData(segments, gen.info, varSimil=0.99, subVariation=0.90)

}

}

\keyword{ utilities }

