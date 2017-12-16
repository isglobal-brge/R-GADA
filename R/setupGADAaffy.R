setupGADAaffy <-
function(file, NumCols, log2ratioCol, MarkerIdCol=1, ChrNameCol=2, ChrPosCol=3, XY=TRUE, sort=FALSE, orderProbes, saveGenInfo=TRUE)
 {  
    if (missing(log2ratioCol))
        stop("Missing log2ratioCol. Please, indicate which column of the file contains the log2ratio. \
              Usually log2ratioCol=4 if exported by the Affymetrix Genotyping Console Tool \
              or log2ratioCol=5 if using affymetrix Power Tools")
    if (missing(NumCols))
        stop("Missing NumCols argument. Please, indicate the number of columns in the file. \
              Usually between 4 and 8 if exported by the Affymetrix Genotyping Console Tool \
              or 8 if using affymetrix Power Tools")
    
    ans<-list() 
    
    xx<-scan(file,what=c("character"),comment.char='#',sep='\t')

    ## Validating headers for the affymetrix files
    headers <- xx[1:NumCols]
    if((headers[MarkerIdCol]!='ProbeSetName')&&(headers[MarkerIdCol]!='ProbeSet'))
      warning("Expecting 'ProbeSetName' as the header of MarkerIdCol in an Affymetrix file")
    if(headers[ChrNameCol]!='Chromosome')
      warning("Expecting 'Chromosome' as the header of ChrNameCol in an Affymetrix file")
    if(headers[ChrPosCol]!='Position')
      warning("Expecting 'Position' as the header of ChrPosCol in an Affymetrix file")
    if(headers[log2ratioCol]!='Log2Ratio')
      warning("Expecting 'Log2Ratio' as the header of  log2ratioCol in an Affymetrix file")   
    
    xx<- xx[-(1:NumCols)] #Remove the header    
    x<-matrix(xx,ncol=NumCols,nrow=length(xx)/NumCols, byrow=TRUE)


    # Removing probes without chromosomes (i.e. Mitocondrial probes)
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

    chrna<-is.na(chr)
    if (sum(chrna)>0)
      {
        warning("Extraneous chromosome names (e.g. 'M') have been removed")
        xx<-x[!chrna,]
        x<-xx   
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
       }

    
    if (!sort)
      {
        ans$log.ratio<-as.numeric(x[,log2ratioCol])
        if(saveGenInfo)
         attr(ans,"gen.info")<-data.frame(probe=x[,MarkerIdCol],chr=chr,pos=as.numeric(x[,ChrPosCol]),stringsAsFactors=FALSE)
        else
         attr(ans,"gen.info")<-NULL
      }

    else  # if sort
      {
        pos<-as.numeric(x[,ChrPosCol])
        if (missing(orderProbes))
          o<-order(chr,pos)
        else
          o<-orderProbes

        ans$log.ratio<-as.numeric(x[o,log2ratioCol])
        if(saveGenInfo)
         attr(ans,"gen.info")<-data.frame(probe=x[o,MarkerIdCol],chr=chr[o],pos=pos[o],stringsAsFactors=FALSE)
        else
         attr(ans,"gen.info")<-NULL
      }   

  attr(ans,"type")<-"Affy" #AffyGTC301
  class(ans)<-"setupGADA"
  ans
 }
