getNamesProbes<-function(x, min.correlation=0.30)
 {
   selec<-x[,6]>=min.correlation
   lab<-x[selec,1]

   if (length(lab)==0)
    {
     warning("No variables selected. Change 'min.correlation' argument")
     ans<-NA
    }

   else
    {
     index<-c("\\.-1","\\.0","\\.1")
     for(i in 1:3)
      lab<-gsub(index[i],"",lab)
 
     cat(" ...", length(lab), "variables selected \n")
     ans<-unique(lab)
    }

   ans
 }
