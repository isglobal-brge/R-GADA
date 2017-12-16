exportToBED<-function(x, num=300, strand="+", col.gain="0,255,255", col.loss="255,0,0", file="BED.txt")
 {  

  if (!inherits(x, "summaryParGADA")) 
   stop("object must be of class 'summaryParGADA'")

  setwd(attr(x,"Info"))

  Samples<-attr(x,"Samples")
  
  dd<-unlist(x,recursive=FALSE) 

  
  resOK<-NULL
  for (k in 1:length(dd))
    {
     resOK<-rbind(resOK,dd[[k]])
    }

  n<-nrow(resOK)
  numOK<-rep(num,n)
  strandOK<-rep(strand,n)
  colorOK<-ifelse(resOK[,6]==-1,col.loss,col.gain)
  chrOK<-paste("chr",resOK[,5],sep="")
  ans<-data.frame(chrOK,resOK[,c(1,2,7)],numOK,strandOK,resOK[,c(1,2)],colorOK)
  
  write.table(ans,file="BED.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  
  cat("File", file, "has been generated at",attr(x,"Info"),"\n") 

  invisible(ans)

 }
