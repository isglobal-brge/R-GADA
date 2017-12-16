`SBL` <-
  function (x, sigma2, aAlpha=0.2, estim.sigma2 = FALSE, maxit = 10000, 
            tol = 1e-08, verbose = FALSE, saveInfo = TRUE) 
{
  if (!inherits(x, "setupGADA")) 
    stop("object must be of class 'setupGADA' ")
  if (missing(sigma2) & !estim.sigma2) 
    stop("parameter 'sigma2' is missing. \n It can be estimated using estim.sigma2=TRUE")
  if (estim.sigma2) 
    sigma2 <- -1
  gen.info <- attr(x, "gen.info")
  
  debug <- ifelse(verbose, 1, 0)

  x$log.ratio[is.na(x$log.ratio)] <- 0  # RPR we finaly use this to deal with the NA?

  if (is.null(gen.info)) {
    out <- list()
    y <- x$log.ratio
    M <- length(y)
    if (sigma2 < 0) {
      sigma2 <- 0.5 * mean(diff(y)^2, trim = 0.01)
      cat("    The estimated sigma2 =", sigma2, "\n")
    }
    out <- SBL.fit(0, y, NULL, NULL, NULL, sigma2, aAlpha, 
                   maxit, tol, debug)
  }
  else {
    chr <- unique(gen.info$chr)
    chr <- chr[!is.na(chr)]

    control <- 1:nrow(gen.info)
    if (sigma2 < 0) {
      sigma2 <- 0
      totalM <- 0
      if (length(chr) > 22) 
        n.chr <- 22
      else n.chr <- length(chr)
      for (i in 1:n.chr) {

# Changed to control missing data, JRG
#        selec <- control[gen.info$chr == chr[i]]

        tt<-gen.info$chr==chr[i]
        tt<-tt[!is.na(tt)]  # RPR:  What is this doing?
        selec <- control[tt]

        y <- x$log.ratio[selec]
        M <- length(y)
        totalM <- totalM + (M - 1)
        sigma2 <- sigma2 + (M - 1) * 0.5 * mean(diff(y)^2, 
                                                trim = 0.01)
      }
      sigma2 <- sigma2/totalM
      cat("    The estimated sigma2 =", sigma2, "\n")
    }
    out <- lapply(1:length(chr), SBL.fit, x = x$log.ratio, 
                  control = control, gen.info = gen.info, chr = chr, 
                  sigma2 = sigma2, aAlpha = aAlpha, maxit = maxit, 
                  tol = tol, debug = debug)
    attr(out, "chr") <- chr  #RPR  I moved this line inside here bc if no gen
  }

  # RPR Potential BUG!!! Do we also account for the removed NA probes in the plots?
  
  attr(out, "data") <- x$log.ratio
  if (saveInfo) {
    attr(out, "gen.info") <- attr(x, "gen.info")
  }
  class(out) <- "SBL"
  out
}
