parBE<-function(x, Samples, T, MinSegLen, verbose=TRUE, ...)
 {

  setwd(x)
 
  if (missing(Samples))
   {
     Samples<-attr(x,"Samples") 

# to control error in parSBL
     temp1<-dir("SBL")
     temp2<-temp1[grep("sbl",temp1)]
     indiv<-gsub("sbl","",temp2) 
     if (length(indiv)!=length(Samples))
      ok<-sort(as.numeric(indiv))
     else
      ok<-NULL
   } 

  if (verbose)
    cat("Retrieving annotation data ...")

  load("SBL/gen.info.Rdata")

  if (verbose)
    cat("done \n")

  if (length(Samples)>2)
   stop(" 'Samples' must be the number of samples or a vector indicating the first and last sample")

  if (length(Samples)==1)
    Samples<-c(1,Samples)
  
  analize.i<-function(i,T, MinSegLen, gen.info, labels, verbose, ok)
    {
      if (verbose)
       cat("   Array #",i,"... ")  
   
      if (i%in%ok | is.null(ok))
       {
        load(paste("SBL/sbl", i, sep = ""))
        attr(step1,"gen.info")<-gen.info
        step2<-BackwardElimination(step1, T=T, MinSegLen=MinSegLen)
        ans<-summary(step2,print=FALSE)
        save(ans,step2,file=paste("SBL/segments",i,sep=""),compress=TRUE)
        ans$sample<-labels[i]
       }
      else
       { 
        ans<-NULL
        ans$sample<-labels[i]
        save(ans, file=paste("SBL/segments",i,sep=""),compress=TRUE)
       } 
      ans     
    }

   if (verbose)
     cat("Backward elimination procedure for",Samples[2]-Samples[1]+1,"samples ... \n")

   res<-plapply(Samples[1]:Samples[2], analize.i, T=T, MinSegLen=MinSegLen, gen.info=gen.info, labels=attr(x,"labels.samples"), verbose=verbose, ok=ok)

   if (verbose)
     cat("Backward elimination procedure for",Samples[2]-Samples[1]+1,"samples ...done \n")

   error<-sum(unlist(lapply(res, function(x) inherits(x, "try-error"))))

   if (error>0)
    {
      cat("WARNING!!! \n")
      cat("  Backward Elimination procedure failed for",sum(error),"samples \n")
      cat("  (type error to see what happened) \n")
      error <<- res
    }

  save(res,file="SBL/allSegments",compress=TRUE)
 }

