\name{setupGADAgeneral}
\alias{setupGADAgeneral}


\title{Imports general array data to gada}

\description{
This function encapsulates information obtained from any platform to be analyzed with gada
}

\usage{
setupGADAgeneral(log.ratio, genotype = NULL, B.allele.freq = NULL, gen.info = NULL, type = "log.ratio", sort = FALSE)
}

\arguments{
  \item{log.ratio}{A vector with log2ratio intensities}
  \item{genotype}{A vector with genotypes. Not required }
  \item{B.allele.freq}{A vector with B-allele frequency. Not required}
  \item{gen.info}{A vector with annotation data. Not required. See details}
  \item{type}{Type of data. By default is 'log.ratio'}
  \item{sort}{Should data be sorted using annotation data given in 'gen.info'}
}

\details{
Annotation data if available can be added through the argument 'gen.info' as a data.frame. The following format is required:

\preformatted{
        probe chr    pos
1  rs12354060   1  10004
2   rs6650104   1 554340
3  rs12184279   1 707348
4  rs12564807   1 724325
5   rs3115860   1 743268
6   rs7515489   1 758845
7  rs17160939   1 773886
8  rs12086311   1 798632
9   rs4475691   1 836671
10 rs28705211   1 890368
...}

}

\value{
An object of class 'setupGADA'
}

\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}

\seealso{ \code{\link{setupGADAIllumina}}, \code{\link{setupGADAaffy}}}

\examples{
\dontrun{
## Simulated data
set.seed(123456)
cn<-rep(c(rep(1,1E5-100),rep(1,100),rep(1,1E5)),4) #Underlying copy number
arrayData<-rnorm(length(cn),mean=(log2(cn)-1),sd=1) #Simulated array

gen.info <- data.frame( probe=paste("id",1:length(cn),sep=""),
            chr=c(rep(1,2E5),rep(2,2E5),rep(3,2E5),rep(4,2E5)),
            pos=rep(1:(length(cn)/4),4)*10)

## setupGADA object with annotation information
dataSim<-setupGADAgeneral(arrayData,gen.info=gen.info)
dataSim
}

}

\keyword{utilities}

