sortPopData<-function(segments,nchr,npop)
{

  ans<-list()
  for(pop in 1:npop)
  {
    rr1<-list() 
    for(chr in 1:nchr) 
    { 
      rr1[[chr]]<-segments[[chr]]
    }
    
    ans[[pop]]<-rr1
   }
ans
} 