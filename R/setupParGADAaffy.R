setupParGADAaffy<-function(folder, files, verbose=TRUE, ...)
 {
 
  if (verbose)
   {
    cat("\n")
    cat("Creating objects of class setupGADA for all input files... \n")    
   }
  
  if(missing(folder))
    folder<-getwd()
  oldpwd<-getwd()
  setwd(folder)

  if (!"rawData"%in%dir())
   stop("rawData folder with the data .txt files cannot be located")  

  if (missing(files))
    files <- dir("rawData",".txt$")
  
  if (!"SBL"%in%dir() )
   dir.create("SBL")

 if (verbose)
  {
  cat("  Applying setupGADAaffy for", length(files) ,"samples ... \n")
  } 


   prepare.i<-function(i, files, ...)
    {
      if (verbose)
       cat("  Importing array: ",files[i],"... ")  

      if (i==1) 
       temp<-setupGADAaffy(paste("rawData/",files[i],sep=""),...)
      else
       temp<-setupGADAaffy(paste("rawData/",files[i],sep=""), saveGenInfo=FALSE, ...)

      save(temp,file=paste("SBL/setupGADA",i,sep=""),compress=TRUE)
     
      if (verbose)
       cat("   Array #",i,"...done \n")  
    }


  res<-plapply(1:length(files),function(i) try(prepare.i(i,files,...), TRUE))

  load("SBL/setupGADA1")
  
  gen.info<-attr(temp,"gen.info");

  save(gen.info,file="SBL/gen.info.Rdata",compress=TRUE)

  attr(temp,"gen.info")<-NULL
  save(temp,file="SBL/setupGADA1",compress=TRUE)

 if (verbose)
  {
   cat("Creating objects of class setupGADA for all input files... done \n")
  }
   error<-sum(unlist(lapply(res, function(x) inherits(x, "try-error"))))

   if (error>0)
    {
      cat("WARNING!!! \n")
      cat("  Creating objects procedure failed for",sum(error),"samples \n")
      cat("  (type error to see what happened) \n")
      error <<- res
    }

 ans<-getwd()
 class(ans)<-"parGADA"
 attr(ans,"type")<-"Affy"
 attr(ans,"labels.samples")<-gsub("sample.","",gsub(".txt","",files)) 
 attr(ans,"Samples")<-length(files)
 ans
}


