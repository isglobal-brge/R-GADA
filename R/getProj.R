
getProj<-function(mat, pop, select=NULL, saveInfo=TRUE, type="acm",verbose=TRUE,...)
{
  
  .GlobalEnv$md<-mat
  .GlobalEnv$pd<-pop

  if (type=="acm")
   {
   
    if (verbose)
    cat("  running acm \n")
   
    if (is.null(select))
     .GlobalEnv$md<-acm.disjonctif(.GlobalEnv$md)
    else
     .GlobalEnv$md<-acm.disjonctif(.GlobalEnv$md)[select]

    .GlobalEnv$dual<-dudi.coa(.GlobalEnv$md,scan=FALSE,...)
   }

  if (type=="pca")
  {
   if (verbose)
    cat("  running pca \n")
	
    dual<-dudi.pca(.GlobalEnv$md,scan=FALSE,...)
  }
  
  if (type=="coa")
   {
   
    if (verbose)
      cat("  running acm \n")
	
	
    if (is.null(select))
     .GlobalEnv$md<-.GlobalEnv$md
    else
     .GlobalEnv$md<-(.GlobalEnv$md)[select]

    .GlobalEnv$dual<-dudi.coa(.GlobalEnv$md,scan=FALSE,...)
   }

  if (verbose)
    cat("  running discrimin \n")
	
  .GlobalEnv$dd<-discrimin(.GlobalEnv$dual,.GlobalEnv$pd,scannf=FALSE,...)
  ans<-project.pop(dd=.GlobalEnv$dd,verbose=verbose)

  if(saveInfo)
   {
    attr(ans,"dudi.dis")<-.GlobalEnv$dd
   }

  return(ans)
}
