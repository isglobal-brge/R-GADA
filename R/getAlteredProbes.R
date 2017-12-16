getAlteredProbes<-function(x, chr, min.perc=0.10, max.number.cnv=100, length.base)
 {

  if (missing(chr))
   chr<-c(1:22)

  if(!inherits(x,"summaryParGADA"))
   stop("object must be of class 'summaryParGADA'")

  if (missing(length.base))
    length.base<-attr(x,"length.base")

  nSamples<-attr(x,"Samples")[2]


  Info<-attr(x,"Info")
  if (is.character(Info))
   {
    load(paste(Info,"/SBL/gen.info.Rdata",sep=""))
   }

  chr.lab<-attr(x, "chr")
  ans<-countAltered.i(chr, x, max.number.cnv, length.base, gen.info, chr.lab)

  cc<-ceiling(nSamples*min.perc) 

  gains<-ans$gains[ans$gains$Freq>=cc,]
  losses<-ans$losses[ans$losses$Freq>=cc,]

  ans<-list(gains=gains, losses=losses)
  ans 
 }
