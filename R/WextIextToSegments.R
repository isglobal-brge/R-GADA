`WextIextToSegments` <-
function(x)
 {
   #Assuming that x is an object with a single chromosome....
   
   Wext<-x$Wext
   Iext<-x$Iext

   K<-length(Wext)-1

   ret<-.C("RcallWextIextToSegments",
	       Wext=as.double(Wext),
               Iext=as.integer(Iext),
	       K=as.integer(K),
	       SegAmp=double(K+1),
	       SegLen=integer(K+1),PACKAGE="gada")

    ans<-data.frame(
	       IniProbe=(ret$Iext[1:(K+1)]+1),
	       EndProbe=ret$Iext[2:(K+2)],
	       LenProbe=ret$SegLen,
	       MeanAmp=ret$SegAmp
                   )
    ans 
}
