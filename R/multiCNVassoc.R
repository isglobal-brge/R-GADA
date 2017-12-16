multiCNVassoc<-function(x,formula,...)
 {
   ans<-plapply(x, function(i) try(assocCNV.i(i,formula,...),TRUE))
   class(ans)<-"multiCNVassoc"
   ans
 }

