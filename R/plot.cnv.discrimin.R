plot.cnv.discrimin<-function(x, all=FALSE, xax = 1, yax = 2, ...)
{
	if (!inherits(x, "cnv.discrimin")) 
        stop("Use only with 'cnv.discrimin' objects")

	.GlobalEnv$pd<-x$pop
	.GlobalEnv$md<-x$mat.f
        
        if (all)
         plot(attr(x$proj,"dudi.dis"))

        else
         {
          xx<-attr(x$proj,"dudi.dis")
          appel <- as.list(xx$call)
          fac <- eval(appel$fac, sys.frame(0))
           s.class(xx$li, fac, xax = xax, yax = yax, sub = "Scores and classes", csub = 2, clab = 1.5)
         }


}
