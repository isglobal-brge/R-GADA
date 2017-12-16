`plotRatio.parGADA` <-
function(x, Sample, segments=FALSE, ...)
{
 setwd(x)
 load("SBL/gen.info.Rdata")

 if (!segments) 
  {
   load(paste("SBL/setupGADA",Sample,sep=""))
   attr(temp, "gen.info")<- gen.info
   plotlogRatio(temp, ...)
  }
 
 else
  {
   load(paste("SBL/segments",Sample,sep=""))
   attr(step2, "gen.info")<- gen.info
   plotRatio(step2, ...)
  }

 title(paste("Sample",Sample), line=-.5)

}

