`print.setupGADA` <-
function(x,...)
 {

  if (!inherits(x, "setupGADA")) 
        stop("object must be of class 'setupGADA' ")

  type<-attr(x,"type")
  gen.info<-attr(x,"gen.info")

  cat("Object of class 'setupGADA' (",type, " data) \n", sep="")
  cat("-------------------------------------------- \n")
  cat("  Number of probes: ", length(x$log.ratio), " (",sum(is.na(x$log.ratio)), " missing values) \n", sep="")
 
  cat("\n")
  cat("  Number of probes by chromosome:")
  if (!is.null(gen.info))
    print(table(gen.info[,2]))
  else
    cat("\n    No genetic information available") 
  

  if (!is.null(x$genotype))
   {
    cat("\n")
    cat("  Genotype frequency: \n")
    tt<-table(x$genotype)
    tt.p<-round(prop.table(tt)*100,1)
  
    out<-NULL
    for (i in 1:length(tt))
     {
      out<-c(out,paste(tt[i]," (",tt.p[i],"%)",sep=""))
     } 
    names(out)<-names(tt)
    print(out,quote=FALSE)
   }


  if (!is.null(x$B.allele.freq))
   {
    cat("\n")
    cat("  Summary of B-allele frequency: \n")
    print(summary(x$B.allele.freq))
   }
   

  cat("\n")
 
 }

