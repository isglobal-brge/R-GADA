rankVariables<-function(x, pop,cnv.blocks=NULL, proj=TRUE, ...)
{
  if (inherits(x,"cnv.discrimin"))
   {
    cp<-x
   }
  else
   {
    cp<-discriminCNV(x, pop,...)
   }

  pr<-sapply(1:length(cp$proj), function(x) cp$proj[[x]][,2])
  psit<-pr<0
  pr[psit]<-0
  
  pr.pop.max<-apply(pr,1,max)
  pr.wh.max<-apply(pr,1,which.max)
  pop.nm<-names(cp$proj)
  pop.nm<-pop.nm[pr.wh.max]
  
  
  mm<-cp$proj[[1]][,1]
  mm<-matrix(as.vector(mm), ncol=1, nrow=length(mm))
  

  if (proj)
   {
    o<-order(-pr.pop.max)
	
	probe<-mm[o]
	
	vm<-as.vector(sapply(probe, function(x) strsplit(x,"\\.")[[1]][1]))
	
    t.cnv.blocks<-data.frame(t(cnv.blocks[,-1]))
	names(t.cnv.blocks)<-as.vector(cnv.blocks[,1])
	
	pr.info<-as.data.frame(sapply(1:length(vm), function(x) t.cnv.blocks[vm[x]]))
	pr.info<-as.data.frame(t(pr.info))
	names(pr.info)<-names(cnv.blocks)[-1]
	rownames(pr.info)<-c()
	
	if(is.null(cnv.blocks))
      ans<-data.frame(probe=probe,correlation=round(pr.pop.max[o],3),population=pop.nm[o])
	else
	  ans<-data.frame(probe=probe,pr.info, correlation=round(pr.pop.max[o],3),population=pop.nm[o])  
   }
  else
   {
    #select variable with axis of MCDA with nrep=1 note: no abs
    o<-order(-apply((attr(cp$proj,"dudi.dis")$va),1,max))
	
	if(missing(cnv.blocks))
       ans<-data.frame(probe=rownames(attr(cp$proj,"dudi.dis")$va)[o], correlation=round(attr(cp$proj,"dudi.dis")$va[o],3),population=pop.nm[o])
	else
	   ans<-data.frame(probe=rownames(attr(cp$proj,"dudi.dis")$va)[o], correlation=round(attr(cp$proj,"dudi.dis")$va[o],3),population=pop.nm[o])   
	   #prob.info non implmented yet for this case!
   }

  if (nrow(ans)>30)
   selec<-30
  else
   selec<-nrow(ans)

  dotchart(rev(ans[1:selec,6]), rev(ans[1:selec,1]), xlab="correlation")
 

 invisible(ans)

}
