setupParGADA <- function(folder, files, MarkerIdCol=1,
                         ChrNameCol=2, ChrPosCol=3,
                         chrs=c(1:22, "X", "Y"),
                         XY=TRUE, verbose=TRUE, 
                         mc.cores=1, ...)
 {
 
  if (!"SBL"%in%dir() )
   dir.create("SBL")
  
  if(missing(folder))
    folder<-getwd()
  oldpwd<-getwd()
  setwd(folder)

  if (!"rawData"%in%dir())
   stop("rawData folder with the data files cannot be located")  

  if (missing(files))
    files <- dir("rawData")

  gen.info <- fread(file.path("rawData", files[1]), 
                    select=c(MarkerIdCol, 
                             ChrNameCol, 
                             ChrPosCol),
                    data.table = FALSE)
  
  colnames(gen.info) <- c("marker", "chr", "position")

  if (XY){
    gen.info$chr[gen.info$chr=="XY"] <- "X"
  }  
  
  gen.info <- gen.info[gen.info$chr%in%chrs,]
  gen.info$chr <- factor(gen.info$chr, 
                         levels= as.character(chrs))
  
  o<-order(match(gen.info$chr, chrs), gen.info$position)
  gen.info<-gen.info[o,]
  attr(gen.info,"sort") <- TRUE
  attr(gen.info,"orderProbe") <- o
  select <- rownames(gen.info)

  save(gen.info, file="SBL/gen.info.Rdata", compress=TRUE)
   

  if (verbose)
   {
    cat("Creating object with annotation data ...done \n")
   }


 if (verbose)
  {
    cat("\n")
    cat("Creating objects of class setupGADA for all input files... \n")
  }


 if (verbose)
  {
  cat("  Applying setupGADAIllumina for", length(files) ,"samples ... \n")
  }



  prepare.i<-function(i, files, select, ...)
    {
      if (verbose)
       cat("  Importing array: ",files[i],"... ")  

      dd <- paste("rawData/",files[i],sep="")

      temp<-setupGADA(dd, orderProbes=select, 
                      saveGenInfo=FALSE, ...)

      save(temp, file=paste("SBL/setupGADA",i,sep="" ), compress=TRUE)

      if (verbose)
       cat("   Array #",i,"...done \n")  
    }


 res <- mclapply(1:length(files), 
              function(i) try(prepare.i(i, files=files, XY=XY, ...), TRUE),
              mc.cores=mc.cores)

 if (verbose)
  {
   cat("  Applying setupGADAIllumina for", length(files) ,"samples ... done \n")
   cat("Creating objects of class setupGADA for all input files... done \n")
  }

   fail<-unlist(lapply(res, function(x) inherits(x, "try-error")))
   error<-sum(fail)

   if (error>0)
    {
      cat("WARNING!!! \n")
      cat("  Creating objects procedure failed for",sum(error),"samples \n")
      cat("  (type error to see what happened) \n")
      cat(paste("    ", error, "files have been removed from the analysis \n"))
      error <<- res
      class(error)<-"error.gada.setup"
    }


 ans <- getwd()
 class(ans) <- "parGADA"
 attr(ans,"labels.samples") <- gsub("sample.","",gsub(".txt","",files))[!fail]
 attr(ans,"Samples") <- length(attr(ans,"labels.samples"))
 ans

}


