findNormalLimits <- function (x) 
{
    setwd(x)

    load("SBL/allSegments")
  
    ff<-function(x,chr=1)
      {
            BaseAmp<-attr(x,"BaseAmp");
            return(BaseAmp[chr])
      }


    intensities<-unlist(lapply(.GlobalEnv$res,FUN=ff))
    intensitiesX<-unlist(lapply(.GlobalEnv$res,FUN=ff,chr="X"))
    intensitiesY<-unlist(lapply(.GlobalEnv$res,FUN=ff,chr="Y"))

    NormIntX<-intensitiesX-intensities;
    SexThreshold<-median(NormIntX);
    EstimatedSex<-(NormIntX<SexThreshold); #True = Male, False = Female

    CopyNumber1BaseAmp<-median(NormIntX[EstimatedSex])-median(NormIntX[!EstimatedSex]);

    limits<-c(CopyNumber1BaseAmp/2,-CopyNumber1BaseAmp/2*(log2(3)-1))


    limits

}
