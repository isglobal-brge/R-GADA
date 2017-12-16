summary.parGADA<-function(object, Samples, threshold, length.base, chr=c(1:22,"X","Y"), ...)
 {  

  x<-object
  setwd(x)
  

  if (missing(length.base))
   { 
     warning("All segments are reported. If you want to filter the minimum and maximum \n  lenght of segments, adjust 'length.base' \n (e.g. length.base=c(500,10e6) in base units)")
     length.base<-c(0,Inf)
   }
   

  if (missing(threshold))
   {
     threshold<-findNormalLimits(x)   
     if(any(is.na(threshold)))
       stop("Normal Limits cannot be estimated. Give 'threshold' argument manually")
   }



  if (missing(Samples))
    Samples<-attr(x,"Samples")

  if (length(Samples)==1)
    Samples<-c(1,Samples)

  load("SBL/allSegments")
  
 
  ff<-function(x,chr,threshold,length.base)
   {
    cond<-x[,5]==chr & x[,6]!=0 & (x[,4]<threshold[1] | x[,4]>threshold[2]) &   x[,2]-x[,1]>=length.base[1] & x[,2]-x[,1]<=length.base[2] & x[,2]-x[,1]>0

    return(x[cond,])
   } 

  ff2<-function(x,chr,threshold,length.base)
   {
    cond<-x[,5]==chr & x[,6]!=0 & x[,2]-x[,1]>0
    xx<-x[cond,]
    cond2<- (xx[,4]<=threshold[1] | xx[,4]>=threshold[2]) &   xx[,2]-xx[,1]>=length.base[1] & xx[,2]-xx[,1]<=length.base[2] 

    return(!cond2)
   } 


  ans<-list()
  no.cnv<-list()
  for (i in 1:length(chr))
   {
   ans[[i]] <-lapply(res,FUN=ff,chr=chr[i],threshold=threshold,length.base=length.base)
   no.cnv[[i]] <-lapply(res,FUN=ff2,chr=chr[i],threshold=threshold,length.base=length.base)
   }

  attr(ans,"no.cnv")<-no.cnv
  attr(ans,"length.base")<-length.base
  attr(ans,"threshold")<-threshold
  attr(ans,"Info")<-x
  attr(ans,"Samples")<-Samples
  attr(ans,"labels.samples")<-labels(x)
  attr(ans,"chr")<-chr
  class(ans)<-"summaryParGADA"
  names(ans)<-paste("chromosome",chr)
  ans
 }

