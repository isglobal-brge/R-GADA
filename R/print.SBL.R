print.SBL <-function(x,...)
 {

  gen.info<-!is.null(attr(x,"gen.info"))

  if (gen.info)
   {
      Info<-collapseInfo(x)
      ans<-data.frame(chromosome=attr(x,"chr"),discontinuities=Info$K,numit=Info$numit,tolerance=Info$tol)    
      cat("Sparse Bayesian Learnig (SBL) algorithm \n")
      cat("sigma2 =",round(x[[1]]$sigma2,4),"\n")
      cat("-------------------------------------------------------------- \n")
      print(ans)   

   }

  else
   {
    cat("Sparse Bayesian Learnig (SBL) algorithm \n")
    cat("SBL coverged after", x$numit ,"iterations with tolerance", x$tol, "\n")
    cat("Total number of discontinuities: ", x$K,"\n") 
    cat("sigma2=",round(x$sigma2,4),"\n")
   }

 }

