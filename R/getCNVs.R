getCNVs <- function(x){
  if(!inherits(x, "summaryParGADA"))
    stop("x must be an object of class summaryParGADA")
  
  chrs <- attr(x, "chr")
  
  ans <- list()
  for (i in chrs)
   ans[[i]] <- rbindlist(x[[paste("chromosome", i)]])
  out <- rbindlist(ans)
  
  out.gr <- GenomicRanges::GRanges(seqnames = paste0("chr", out$chromosome),
                    ranges = IRanges::IRanges(
                      start = out$IniProbe,
                      end = out$EndProbe),
                    LenProbe = out$LenProbe,
                    MeanAmp = round(out$MeanAmp, 2),
                    State = out$State, 
                    sample = out$sample)
  out.gr
}
