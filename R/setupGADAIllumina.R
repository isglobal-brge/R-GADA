setupGADAIllumina <-
function(file, NumCols, log2ratioCol, MarkerIdCol=1, ChrNameCol=2, ChrPosCol=3, XY=TRUE, sort=TRUE, orderProbes, saveGenInfo=TRUE)
 {
 
    if (missing(log2ratioCol))
        stop("Missing log2ratioCol. Please, indicate which column of the file contains the log2ratio")
    if (missing(NumCols))
        stop("Missing NumCols argument. Please, indicate the number of columns in the file")
    
    ans<-list() 
    xx<-scan(file, skip=1, what=c("character"), sep="\t")


    ## Validating headers for the Illumina files
    headers <- scan(file, nline=1, what=c("character"), quiet=TRUE, sep="\t")

    if(headers[MarkerIdCol]!='Name')
      warning("Expecting 'Name' as the header of MarkerIdCol in an Illumina file")
    if(headers[ChrNameCol]!='Chr')
      warning("Expecting 'Chr' as the header of ChrNameCol in an Illumina file")
    if(!headers[ChrPosCol]%in%c('position', 'Position'))
      warning("Expecting 'position' as the header of ChrPosCol in an Illumina file")
    if(NROW(grep("Log.R.Ratio",headers[log2ratioCol]))!=1)
      warning("Expecting 'Log R Ratio' as the header of  log2ratioCol in an Illumina file")   

    x<-matrix(xx, ncol=NumCols, nrow=length(xx)/NumCols, byrow=TRUE)
  
    if (!sort)
     {
      x[,ChrNameCol][x[,ChrNameCol]=="XY"]<-"X"
      if (XY)
       {
        chr<-factor(x[,ChrNameCol],levels=c(as.character(1:22),"X","Y"))
       }
      else
       {
        x[,ChrNameCol][x[,ChrNameCol]=="23"]<-"X"
        x[,ChrNameCol][x[,ChrNameCol]=="24"]<-"Y"
        chr<-factor(x[,ChrNameCol],levels=c(as.character(1:22),"X","Y"))
       }

      temp<-data.frame(probe=x[,MarkerIdCol],chr=chr,pos=as.numeric(x[,ChrPosCol]),stringsAsFactors=FALSE)

# Mitocondrial (or missing)
      mito<-is.na(temp$chr) | is.na(temp$pos)
      temp2<-temp[!mito,] 

      ans$log.ratio<-as.numeric(x[!mito,log2ratioCol])
      attr(ans,"gen.info")<-temp2
     }

    else  # if sort
     {
      x[,ChrNameCol][x[,ChrNameCol]=="XY"]<-"X"
      if (XY)
       {
        chr<-factor(x[,ChrNameCol],levels=c(as.character(1:22),"X","Y"))
       }
      else
       {
        x[,ChrNameCol][x[,ChrNameCol]=="23"]<-"X"
        x[,ChrNameCol][x[,ChrNameCol]=="24"]<-"Y"
        chr<-factor(x[,ChrNameCol],levels=c(as.character(1:22),"X","Y"))
       }

      pos<-as.numeric(x[, ChrPosCol])
      if (missing(orderProbes))
       o<-order(chr,pos)
      else
       o<-orderProbes

      temp<-data.frame(probe=x[o,MarkerIdCol],chr=chr[o],pos=pos[o],stringsAsFactors=FALSE)

# Mitocondrial
      mito<-is.na(temp$chr) | is.na(temp$pos)
      temp2<-temp[!mito,] 

      aux<-as.numeric(x[o,log2ratioCol]) 
      ans$log.ratio<-aux[!mito]
      attr(ans,"gen.info")<-temp2
     }

  attr(ans,"type")<-"Illumina"
  if (!saveGenInfo)
   attr(ans, "gen.info")<-TRUE

  class(ans)<-"setupGADA"
  ans
 }

