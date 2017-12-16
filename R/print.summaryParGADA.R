print.summaryParGADA<-function(x, ...)
 {
  if (!inherits(x, "summaryParGADA"))
   stop("object must be of class 'summaryParGADA' ")

  length.base<-attr(x,"length.base")
  threshold<-round(attr(x,"threshold"),2)
  no.cnv<-attr(x,"no.cnv")

  nInd<-attr(x,"Samples")[2] 

  out<-list()
  CNVs<-list()
  tamanys<-NULL

  chr<-attr(x, "chr")
  n.chr<-length(chr)
  if (n.chr>22)
   n.chr<-22

  for (i in 1:n.chr)
  # No cromosoma X Y -> Luego preparar por cromosoma
  {
   tt<-x[[i]]
   gains <- 0
   losses <- 0
   n.cnv <- 0
   for (j in 1:nInd)
    {
     tt2<-tt[[j]]

     tamanys<-c(tamanys,tt2[,2]-tt2[,1])

     n.cnv.j<-nrow(tt2)
     n.cnv<-n.cnv+n.cnv.j

     gains<-gains+sum(tt2$State==1)
     losses<-losses+sum(tt2$State==-1)
    }

   CNVs[[i]]<-n.cnv

   out[[i]]<-c(n.cnv,gains,losses)

   n.cnvT <- sum(unlist(lapply(out, function(x) x[1])))
   gainsT <- sum(unlist(lapply(out, function(x) x[2])))
   lossesT <- sum(unlist(lapply(out, function(x) x[3])))
   per.gainsT<-round((gainsT/(n.cnvT))*100,1)
   per.lossesT<-round((lossesT/(n.cnvT))*100,1)

   no.cnvT<-sum(unlist(no.cnv))

    res<-list(no.cnvT,c(n.cnvT,gainsT,per.gainsT,lossesT,per.lossesT),summary(tamanys),CNVs,out)
  }
  
  cat("\n")
  cat("------------------------------------------- \n")
  cat("Summary results for", nInd, "individuals \n")
  cat("------------------------------------------- \n")

  cat("NOTE: ", no.cnvT, " segments with length not in the range ", length.base[1], "-", length.base[2], " bases", " and with mean log2ratio in the range (", threshold[1], ",", threshold[2], ") have been discarded \n", sep="")



  total<-data.frame(t(res[[2]]))
  names(total)<-c("# segments","Gains","%","Losses","%")
  dimnames(total)[[1]]<-""

  cat("\n")
  cat("Number of Total Segments:","\n")
  print(total)

  cat("\n")
  cat("Summary of length of segments: ","\n",sep="")
  print(res[[3]])


  ii<-t(data.frame(res[[5]]))
  dimnames(ii)[[1]]<-paste("Chromosome",1:n.chr)
  dimnames(ii)[[2]]<-c("segments", "Gains", "Losses")

  cat("\n")
  cat("Number of Total Segments by chromosome:","\n")
  print(ii)

  cat("\n")


  invisible(res)
 }
