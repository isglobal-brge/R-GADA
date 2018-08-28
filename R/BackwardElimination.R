`BackwardElimination` <-
function(x, T, MinSegLen, saveInfo=TRUE) 
 {

  if (!inherits(x, "SBL") && !inherits(x, "BackwardElimination")) 
   stop("object must be of class 'SBL' or 'BackwardElimination' ")

  if (T<=0) 
   stop("T must be positive")


  gen.info<-attr(x,"gen.info")

  if (is.null(gen.info))
   {
    Wext<-x$Wext
    Iext<-x$Iext 
    sigma2<-x$sigma2

    K<-length(Wext)-1

    out<-list()
    ret<-.C("RcallBEwTandMinLen",
	       Wext=as.double(Wext),
               Iext=as.integer(Iext),
	       K=as.integer(K),         
               sigma2=as.double(sigma2),
               T=as.double(T),
	       MinSegLen=as.integer(MinSegLen),PACKAGE="gada")    
    K<-ret$K
    out$K<-K
    out$T<-ret$T
    out$MinSegLen<-ret$MinSegLen
    out$Iext <- ret$Iext[1:(K+2)];
    out$Wext <- ret$Wext[1:(K+1)];
    
    out$sigma2<-ret$sigma2
   }
  else
   { 
    chr <- attr(x, "chr")
    out<- lapply(1:length(chr), BackwardElimination.fit, x=x, chr=chr, T=T, MinSegLen=MinSegLen)
    attr(out, "chr") <- chr  
   }  

   attr(out,"data")<-attr(x,"data")
   if (saveInfo)
    attr(out,"gen.info")<-attr(x,"gen.info")

   class(out)<-"BackwardElimination"
   out
}

