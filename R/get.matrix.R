get.matrix<-function(crom.data, probe.pos, verbose,size.min=500,size.max=1.e6){
  #names
   rown<-as.character(1:length(crom.data))
   coln<-names(probe.pos)

  #matrix
   Mat<-matrix(0,length(crom.data),length(probe.pos), dimnames=list(rown,coln))

  #loop over subjects
   for(sub in 1:length(crom.data)){
     if (verbose)
     cat(c("   adding subject:", sub, "\n"))

     numprob<-dim(crom.data[[sub]])[1]
     if (numprob!=0)
     {  
      #loop over probes
       for(prob in 1:numprob)
         {
          pos<-c(crom.data[[sub]][prob,1],crom.data[[sub]][prob,2])
		  
          #get the size of probes
		  size.probe<-crom.data[[sub]][prob,2]-crom.data[[sub]][prob,1]
          
		  #get index and state
          g.ind<-get.index(pos,probe.pos)

          state<-crom.data[[sub]]$State[prob]
 
# POS<<-pos
# PROBE<<-probe.pos 
 
                   #loop over indexes within each probe
          for(ind in g.ind[1]:g.ind[2])
           {
              if (size.probe>=size.min & size.probe<= size.max)
                 Mat[sub,ind]<-state     
           }
        }
      }
 }
 cat("  matrix done \n")
 Mat
}
