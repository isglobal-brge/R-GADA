findNormalLimits <- function (x) 
{
  
  local_env <- new.env()
  load(file.path(x, "SBL/allSegments"), envir = local_env)
  
  ff<-function(x, chr=1){
    BaseAmp <- attr(x,"BaseAmp")
    return(BaseAmp[chr])
  }
  
  intensities <- unlist(lapply(get("res", envir = local_env), 
                               FUN=ff))
  intensitiesX <- unlist(lapply(get("res", envir = local_env),
                                FUN=ff, chr="X"))
  intensitiesY <- unlist(lapply(get("res", envir = local_env), 
                                FUN=ff, chr="Y"))
  
  NormIntX <- intensitiesX - intensities
  SexThreshold <- median(NormIntX, na.rm=TRUE)
  EstimatedSex<-(NormIntX<SexThreshold) #True = Male, False = Female
  
  CopyNumber1BaseAmp <- median(NormIntX[EstimatedSex], na.rm=TRUE)-
    median(NormIntX[!EstimatedSex], na.rm=TRUE)
  
  limits<-c(CopyNumber1BaseAmp/2, -CopyNumber1BaseAmp/2*(log2(3)-1))
  limits
  
}
