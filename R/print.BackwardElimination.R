print.BackwardElimination <-function(x,...)
 {

  gen.info<-!is.null(attr(x,"gen.info"))

  if (gen.info)
   {
      Info<-collapseInfo(x)
      ans<-data.frame(chromosome=attr(x,"chr"),discontinuities=Info$K)    
      cat("Sparse Bayesian Learning (SBL) algorithm \n")
      cat("SBL and Backward Elimination with T=",x[[1]]$T," and minimun length size=",x[[1]]$MinSegLen, "\n",sep="")
      cat("sigma2 =",round(x[[1]]$sigma2,4),"\n")
      cat("-------------------------------------------------------------- \n")
      print(ans)   

   }

  else
   {
    cat("Sparse Bayesian Learning (SBL) algorithm \n")
    cat("SBL and Backward Elimination with T=",x$T," and minimun length size=",x$MinSegLen, "\n",sep="")
    cat("SBL coverged after", x$numit ,"iterations with tolerance", x$tol, "\n")

    cat("Total number of discontinuities: ", x$K,"\n") 
  
    cat("sigma2=",round(x$sigma2,4),"\n")
   }

 }

