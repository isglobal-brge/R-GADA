splitDataBeadStudio<-function(file, Samples, NumCols, verbose=TRUE)
 {
 
  if (missing(Samples))
   stop("Argument 'Samples' must be given")

  if (NumCols<=3)
    stop("'NumCols must be greater than 3 since it includes three columns for annotation data")

# First three columns must have annotation data
  NumCols<-NumCols-3

  if (verbose)
   {
    cat("\n")
    cat("Spliting data from BeadStudio ... \n")
   }

  if (!"rawData"%in%dir())
   dir.create("rawData")
  if (!"SBL"%in%dir() )
   dir.create("SBL")
  
  splitData<-system.file("exec/splitData.pl",package="gada")

  if (verbose)
   {
    cat("\n")
    cat("Obtaining Ratio Intensity files ... \n")
   }

  Samples<-Samples+1

  pp<-paste("perl ",splitData,file,"split", NumCols, "-heading 3 -tab  -out rawData/sample --name_by_header")

  system(pp)

  if (verbose)
   {
    cat("Obtaining Ratio Intensity files ... done \n")
   }

}


