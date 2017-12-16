\name{getChromosomeDat}
\Rdversion{1.1}
\alias{getChromosomeDat}

\title{Cromosome CNV Data}

\description{
Obtains a matix of CNV segments and subject samples, 
from segmentation output
}

\usage{
getChromosomeDat(object, ch, gg, verbose = TRUE)
}

\arguments{
  \item{object}{concatenated list containing segmentation outcome for each subject, chromosome and population. }
  \item{ch}{chromosome number.}

  \item{gg}{genetic info: \code{data.frame} with probe names, chromosome number and genomic positions of each CNV.}

  \item{verbose}{whether to print detailed progress on R prompt.}
}

\details{
For \code{object}, concatenation by subject and chromosome, within a single population,
is readily provided by the function \code{summary.parGADA} on an object \code{parGADA} 
that has been previously segemnted. Population concatenation is the result of   
\code{c(summary{parGADA1}), summary{parGADA2}), ...} across groups.
}


\value{
returns a \code{matrix} with segment calling (-1:loss,0:nochange,1:gain) for each probe across the whole subject sample. 
}


\seealso{
\code{\link{reduceMatrix}}
}

\examples{

\dontrun{

# 'HapMap270Cnv.RData' file can be downloaded from ...

load("Hapmap270Cnv.RData")

ch<-6
mat<-get.chromosome.dat(segments,ch,gen.info)
dim(mat)
mat[1:5,1:5]
}

}

\keyword{ utilities }

