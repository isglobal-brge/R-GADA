\name{reduceMatrix}
\Rdversion{1.1}
\alias{reduceMatrix}
\title{
Reduction of matrix of CNV segments
}
\description{
Reduces the columns of matrix of CNV segments by gathering similar neigboring 
colums and discarding those whose values do not vary much across rows.  
}
\usage{
reduceMatrix(x, gen.info, chr, varSimil = 0.99, subVariation = 0.99, inc = 0)
}
\arguments{
  \item{x}{
Matrix of CNV segements obtanied by getChromosomeDat
}
  \item{gen.info}{
genetic info: genomic positions of each probe. 
}
  \item{chr}{
Chromosome number.
}
  \item{varSimil}{
\code{varSimil=0.99} produces blocks of neighboring CNVs which do differ 
in more than 99\% across the subject sample. 
}
  \item{subVariation}{
\code{subVariation = 0.99} gets rid of CNV blocks that are constant across the 
the whole subject sample, except for 1\% of the subjects.
}
  \item{inc}{
scalar to change the coding of losses, gains or no-changes in CNV calling
}
}
\details{
Similarity between variables is taken respect the first variable of the block. 
Therefore, the extreems of all blocks do not differ in more than \code{varSim \%} 
of their entries. The value of the block across subjects is that of the first CNV
in the block. 
}
\value{
\item{reduceMatrix }{ produces a \code{data.frame} with colums \code{as.factors}
representing each CNV block}
\item{attr(reduceMatrix,"cnv.blocks"}{lists the genetic information of each CNV block}
}

\author{
Alejandro Caceres}

\seealso{
\code{\link{getChromosomeDat}},\code{\link{getReducedData}},\code{\link{matrixPlot}}
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

