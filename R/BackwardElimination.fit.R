BackwardElimination.fit<-function(i,x,chr,T,MinSegLen)
  {
    ans<-list()
    
    Wext<-x[[i]]$Wext
    Iext<-x[[i]]$Iext 
    sigma2<-x[[i]]$sigma2
    K<-length(Wext)-1
    
    ret<-.C("RcallBEwTandMinLen",
            Wext=as.double(Wext),
            Iext=as.integer(Iext),
            K=as.integer(K),         
            sigma2=as.double(sigma2),
            T=as.double(T),
            MinSegLen=as.integer(MinSegLen),PACKAGE="gada")
    
    K<-ret$K
    ans$K<-K
    ans$T<-ret$T
    ans$MinSegLen<-ret$MinSegLen
    ans$Iext <- ret$Iext[1:(K+2)];
    ans$Wext <- ret$Wext[1:(K+1)];
    ans$sigma2<-ret$sigma2
    
    class(ans)<-"BackwardElimination"
    attr(ans, "chr") <- chr[i]
    ans 
  }
