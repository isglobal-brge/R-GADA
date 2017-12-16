setupParGADAIllumina<-function(folder, files, XY=TRUE, verbose=TRUE, sort=TRUE, ...)
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

  splitData<-system.file("exec/splitData.pl",package="gada")

  if (verbose)
   {
    cat("\n")
    cat("Creating object with annotation data ... \n")
   }

  pp<-paste("perl ",splitData, paste("rawData/", files[1], sep=""),"-start_split 1 -end_split 1 split 3 -tab -out SBL/gen.info ")
  system(pp)

  file.rename("SBL/gen.info.n1","SBL/genomicInfo")


  if (sort)
   {
    gg<-scan("SBL/genomicInfo",skip=1,what="character")
    gg2<-matrix(gg,ncol=3,nrow=length(gg)/3,byrow=TRUE)
    gg2[,2][gg2[,2]=="XY"]<-"X"

    if (!XY)
     {
      gg2[,2][gg2[,2]=="23"]<-"X"
      gg2[,2][gg2[,2]=="24"]<-"X"
     }

    temp<-data.frame(probe=gg2[,1], chr=factor(gg2[,2],levels=c(1:22,"X","Y")), pos=as.numeric(gg2[,3]),stringsAsFactors=FALSE) 

# mitocondrial?
    mito<-is.na(temp$chr) 
    gen.info<-temp[!mito,]

    attr(gen.info,"sort")<-TRUE
    o<-order(gen.info$chr,gen.info$pos)
    gen.info<-gen.info[o,]
    attr(gen.info,"orderProbe")<-o
   }

  else
   {
    gg<-scan("SBL/genomicInfo",skip=1,what="character")
    gg2<-matrix(gg,ncol=3,nrow=length(gg)/3,byrow=TRUE)
    gg2[,2][gg2[,2]=="XY"]<-"X"

    if (!XY)
     {
      gg2[,2][gg2[,2]=="23"]<-"X"
      gg2[,2][gg2[,2]=="24"]<-"X"
     }

    temp<-data.frame(probe=gg2[,1], chr=factor(gg2[,2],levels=c(1:22,"X","Y")), pos=as.numeric(gg2[,3]),stringsAsFactors=FALSE) 

# mitocondrial?
    mito<-is.na(temp$chr) 
    gen.info<-temp[!mito,]
    attr(gen.info,"sort")<-FALSE
   } 
  
  save(gen.info,file="SBL/gen.info.Rdata",compress=TRUE)
   

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



  prepare.i<-function(i, files, XY=XY, ...)
    {
      if (verbose)
       cat("  Importing array: ",files[i],"... ")  

      dd<-paste("rawData/",files[i],sep="")
      temp<-setupGADAIllumina(dd, XY=XY, saveGenInfo=FALSE, ...)

      save(temp, file=paste("SBL/setupGADA",i,sep="" ), compress=TRUE)

      if (verbose)
       cat("   Array #",i,"...done \n")  
    }


 res<-plapply(1:length(files), function(i) try(prepare.i(i, files=files, XY=XY, ...), TRUE))

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


 ans<-getwd()
 class(ans)<-"parGADA"
 attr(ans,"type")<-"Illumina"
 attr(ans,"labels.samples")<-gsub("sample.","",gsub(".txt","",files))[!fail]
 attr(ans,"Samples")<-length(attr(ans,"labels.samples"))
 ans

}


