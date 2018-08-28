setupGADA <-
  function(file, MarkerIdCol=1, ChrNameCol=2, ChrPosCol=3,
           log2ratioCol, BAFcol, chrs = c(1:22, "X", "Y"),
           XY=TRUE, sort = TRUE, orderProbes, gen.info,
           saveGenInfo = TRUE) 
{
    if (missing(log2ratioCol)) 
      stop("Missing log2ratioCol. Please, indicate which column of the file contains the LRR")
    if (missing(BAFcol)) 
      stop("Missing BAFcol. Please, indicate which column of the file contains the BAF")
   
    ans <- list()
    
    x <- fread(file, data.table = FALSE)
    
    if (XY){
      x[,ChrNameCol][x[,ChrNameCol]=="XY"] <- "X"
    }
    
    x <- x[x[,ChrNameCol]%in%chrs,]
    
    if (!missing(orderProbes)) {
      x <- x[orderProbes,]
    }  
    else {
      if (sort)
        x <- x[order(match(x[,ChrNameCol], chrs), x[, ChrPosCol]), ]
    }
    
    ans$log.ratio <- x[,log2ratioCol]
    nas <- is.na(ans$log.ratio)
    if (mean(nas)>0.05)
      stop("There are more than 5% of missing data in LRR")
    else
      ans$log.ratio[nas] <- 0
    
    ans$B.allele.freq <- x[, BAFcol]
    
    if (saveGenInfo) {
      if (missing(gen.info)){
        gen.info <- x[,c(MarkerIdCol, ChrNameCol, ChrPosCol)]
        gen.info[,2] <- factor(gen.info[,2], 
                               levels= as.character(chrs))
        names(gen.info) <- c("name", "chr", "position")
        rownames(gen.info) <- 1:nrow(gen.info) # important!!! 08/2018
        attr(ans, "gen.info") <- gen.info
      }
      else {  
        rownames(gen.info) <- 1:nrow(gen.info) # important!!! 08/2018
        attr(ans, "gen.info") <- gen.info
      } 
    }    
    else {
     attr(ans, "gen.info") <- NULL  
    } 
    
    class(ans)<- "setupGADA"
    ans
}
