getReducedData<-function(x, gg, varSimil=0.99, subVariation=0.99)
 {
  cat("Warning !!!: Only autosomes are analyzed \n")
  mat<-NULL
  for (i in 1:22)
   {
    temp<-getChromosomeDat(x, i, gg, verbose=FALSE)
    mat.i<-reduceMatrix(temp, gg, i, varSimil, subVariation)
    if (is.null(mat))
	{
     mat<-mat.i
     cnv.blocks<-attr(mat.i,"cnv.blocks")
	 }
    else
     {
      if (!is.null(dim(mat.i)))
       mat<-cbind(mat, mat.i)
       cnv.blocks<-rbind(cnv.blocks,attr(mat.i,"cnv.blocks"))
     }
   }
   attr(mat,"cnv.blocks")<-cnv.blocks
   mat
 }


