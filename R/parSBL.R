parSBL<-function(x, Samples, estim.sigma2, aAlpha=0.2, 
                 verbose=TRUE, mc.cores=1, ...)
 {

  setwd(x)

  if (verbose)
   cat("Creating SBL directory ...")

  if (!"SBL"%in%dir())
   system("mkdir SBL")

  if (verbose)
    cat("done \n")


  if (missing(Samples))
    Samples<-attr(x,"Samples")


  if (length(Samples)>2)
   stop(" 'Samples' must be the number of samples or a vector indicating the first and last sample")

  if (length(Samples)==1)
   Samples<-c(1,Samples)

  if (verbose)
    cat("Retrieving annotation data ...")

  load("SBL/gen.info.Rdata")

  if (verbose)
    cat("done \n")
 


  analize.i<-function(i,estim.sigma2, aAlpha, gen.info, verbose) 
    {
      if (verbose)
       cat("   Array #",i,"...") 

      load(paste("SBL/setupGADA",i,sep=""))
      attr(temp,"gen.info")<-gen.info  
      step1<-SBL(temp, estim.sigma2=estim.sigma2, aAlpha=aAlpha, saveInfo=FALSE)
      save(step1,file=paste("SBL/sbl",i,sep="" ),compress=TRUE)

      if (verbose)
       cat("   Array #",i,"...done \n")       
    }

   if (verbose)
     cat("Segmentation procedure for",Samples[2]-Samples[1]+1,"samples ... \n")

   res <- mclapply(Samples[1]:Samples[2],
                  function(i) try(analize.i(i, estim.sigma2=estim.sigma2, aAlpha=aAlpha, gen.info=gen.info, verbose=verbose), TRUE),
                  mc.cores=mc.cores)


   if (verbose)
     cat("Segmentation procedure for",Samples[2]-Samples[1]+1,"samples ...done \n")


   fail<-unlist(lapply(res, function(x) inherits(x, "try-error")))
   error<-sum(fail)

   if (error>0)
    {
      cat("WARNING!!! \n")
      cat("  Segmentation procedure failed for",sum(error),"samples \n")
      cat("  (type error to see what happened) \n")
      cat(paste("    ", error, "files have been removed from the analysis \n"))
      error <<- res
      class(error)<-"error.gada.sbl"
    }
 }
