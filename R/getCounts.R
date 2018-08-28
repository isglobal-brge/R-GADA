getCounts <- function(x, group, id="sample"){
  
  if(!inherits(x, "GRanges"))
    stop("'x' should be a 'GRanges object created using 'getCNVs' function")
  
  agg <- aggregate(rep(1, length(x)), 
                   by=list(group=mcols(x)[, group], 
                           sample=mcols(x)[, id]), 
                   FUN=sum)
  n <- table(agg$group)
  
  gg <- mcols(x)[, group]
  groups <- unique(gg)
  
  if(length(groups) > 2)
    stop("grouping variable must have only 2 levels")
  
  x.g1 <- x[gg==groups[1]]
  ss1 <- unique(mcols(x.g1)[, id])
  x.g1.list <- GRangesList(lapply(ss1, function(x) 
    subset(x.g1, get(id)==x)))
  counts.g1 <- countOverlaps(x, x.g1.list)
  
  x.g2 <- x[gg==groups[2]]
  ss2 <- unique(mcols(x.g2)[, id])
  x.g2.list <- GRangesList(lapply(ss2, function(x) 
    subset(x.g2, get(id)==x)))
  counts.g2 <- countOverlaps(x, x.g2.list)
  
  counts <- cbind(counts.g1, counts.g2)
  rownames(counts) <- c(paste(seqnames(x.g1),ranges(x.g1), sep=":"),
                        paste(seqnames(x.g2),ranges(x.g2), sep=":"))
  colnames(counts) <- groups
  ans <- list(counts=counts, n=n)
  ans
}
