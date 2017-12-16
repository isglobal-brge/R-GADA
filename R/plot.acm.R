plot.acm<-function(x, pop, var=TRUE, which.axes=c(1,2), pos="topleft", mycol=c("red","blue","green","yellow"),td=FALSE, ll=levels(pop), pnt=0.2,...)
{
 if (!inherits(x, "acm")) 
         stop("Use only with 'acm' objects")
 
 if (missing(pop))
   pop<-factor(rep(1,length(x$li)))
 else
   if (!is.factor(pop))
     stop(" 'pop' must be a factor variable") 
 
 d<-x$li
 c<-x$co
   
 nPop<-length(ll) 

if (td==FALSE)
{
 if (max(which.axes)>NCOL(d))
  stop("Try again 'dudi.acm' and increase the argument 'nf' ")
 
 if(var)
 {
   plot(rbind(as.matrix(d)[,which.axes], as.matrix(c)[,which.axes]), type="n", ...)
   for(i in 1:nPop)
    {
     sel<-pop==ll[i]
     points(d[sel,which.axes], pch= 19, col=mycol[i])
    }
    points(c[,which.axes], pch=21, cex=0.7)
  } 
  else 
 { 
   plot(rbind(as.matrix(d)[,which.axes]), type="n", ...) 
   for(i in 1:nPop)
    {
     sel<-pop==ll[i]
     points(d[sel,which.axes], pch= 19, cex=pnt,col=mycol[i])
    }
  }
 
 abline(h=0, v=0, lty=2)
 legend(pos, ll, pch=19, col=mycol[1:nPop])

}
else
{
  if (dim(d)[2]!=3)
  stop("Run'dudi.acm' for nf=3 ")
  
  vc<-vector(mode = "logical", length = length(pop))
  cl<-1:length(pop)
  
  for(i in 1:nPop)
  {
    sel<-pop==ll[i]  
    vc<-vc | sel
	cl[sel]<-mycol[i]
  }
  
  plot3d(d[vc,],col=cl[vc],size=.8,type="s")
	
}


} 
