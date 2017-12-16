setupGADAgeneral <-
  function (log.ratio, genotype = NULL, B.allele.freq = NULL, gen.info = NULL, 
    type = "log.ratio", sort = FALSE) 
{
    ans <- list()
    ans$log.ratio <- log.ratio
    ans$genotype <- genotype
    ans$B.allele.freq <- ans$B.allele.freq
    attr(ans, "type") <- type
    attr(ans, "gen.info") <- gen.info
    class(ans)<- "setupGADA";
ans
}
