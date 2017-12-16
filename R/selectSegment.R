selectSegment<-function(x,o)
 {
   oo<-c(1:length(o))[o]
   tt1<-x$IniProbe<=oo[1] 
   tt2<-x$EndProbe>=oo[length(oo)] 

   xx<-x[sum(tt1):(sum(!tt2)+1),] 
   if (oo[1]>xx$IniProbe[1])
    {
     dif<-xx$EndProbe[1]-oo[1]
     xx$LenProbe[1]<-dif+1
    }
      
   if (oo[length(oo)]<xx$EndProbe[nrow(xx)])
    {
     dif<-oo[length(oo)]-xx$IniProbe[nrow(xx)]
     xx$LenProbe[nrow(xx)]<-dif+1
    }
  xx 
 }
