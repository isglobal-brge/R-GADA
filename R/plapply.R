plapply<-function(X, FUN, ...){
if (exists("parLapply") & exists(".GlobalEnv$cl"))
  {o<-parLapply(.GlobalEnv$cl,X, FUN,...)}
else
  {o<-lapply(X, FUN,...)}
o
}
