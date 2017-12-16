SBL.fit<-function(i,x,control,gen.info,chr,sigma2,aAlpha,maxit,tol,debug)
  {              
    ans<-list()
    if (i==0) #  analysis using whole data   
      {
        y<-x
        M <- length(y)
      }
    else      # analysis chromosome by chromosome 
      {  
        ## Changed to control missing data, JRG
        ##                   selec <- control[gen.info$chr == chr[i]]
        
        tt<-gen.info$chr==chr[i]
        tt<-tt[!is.na(tt)]
        selec <- control[tt]
        y <- x[selec]
        M <- length(y) 
      } 

    ret <- .C("Rcall_Single_SBL_PWC_norm", y = as.double(y), 
              Iext = as.integer(0:M), alpha = double(M - 1), 
              Wext = double(M), Wsigma2 = double(M - 1), K = as.integer(M - 1),
              sigma2 = as.double(sigma2), aAlpha = as.double(aAlpha), 
              bAlpha = as.double(0), maxAlpha = as.double(1e+08), 
              isRemovingOK = as.integer(1), numit = as.integer(maxit), 
              tol = as.double(tol), debug = as.integer(debug), 
              PACKAGE = "gada")

    K <- ret$K
    ans$K <- K
    ans$T <- 0
    ans$MinSegLen <- 0
    ans$Iext <- ret$Iext[1:(K + 2)]
    ans$Wext <- ret$Wext[1:(K + 1)]
    ans$sigma2 <- ret$sigma2
    ans$numit <- ret$numit
    ans$tol <- ret$tol
    class(ans) <- "SBL"
    attr(ans, "chr") <- chr[i]
    if ((ret$numit >= maxit) || (ret$tol >= tol))
      warning(paste("SBL algorithm did not converge after ", maxit, "iterations and change", ret$tol))
    ans

  }

