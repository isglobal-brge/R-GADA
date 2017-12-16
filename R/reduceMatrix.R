reduceMatrix<-function(x, gen.info, chr, varSimil=0.99, subVariation=0.99, inc=0)
{
	cat("  reducing matrix... \n")
	
	if(names(gen.info)[2]!="chr")
		warning("getting chromosome number from second column of gen.info")
	
	gg<-gen.info[gen.info[,2]==chr,]
	chr.vect<-rep(chr,dim(x)[2])
	nux<-x+inc		
	
	#gather variables into groups where they do not differ more than "varSimil%" across subjects
	
	#gather variables into groups where they do not differ more than "varSimil%" across subjects

	sub<-1
	s1<-1
	for(c in 1:(dim(nux)[2]-1)) 
	{ 
		if(chr.vect[s1]!=chr.vect[c+1] | sum(nux[,s1]==nux[,c+1])<floor(dim(nux)[1]*varSimil) )  
		{
			s1<-c+1
			sub<-c(sub,s1)
		}
		
	}
	ss<-sub
	
	block.sup<- c((ss-1)[-1],dim(nux)[2])
	block.inf<- ss

	cnv.blocks<-data.frame(num.pr=(block.sup-block.inf+1), pos.inf=gg[block.inf,3],pos.sup=gg[block.sup,3],chr=gg[block.inf,2])
	
	#take probes that have at least less than "subVariation%" of subjects in their most frequent population	
	x.f<-data.frame(apply(nux[,ss],2,as.factor))
	
	num.fac<-lapply(1:length(x.f), function(y) table(x.f[,y]))
	max.num.fac<-unlist(lapply(num.fac,function(y) max(unlist(y))))
	
	varied.probes<-max.num.fac<=floor(dim(x.f)[1]*subVariation)

	cnv.blocks<-cnv.blocks[varied.probes,]

	numPr<-sum(as.numeric(varied.probes))

	if (numPr==0)
		x.f<-c()
	else
     {
	if (numPr==1)
		x.f<-data.frame(array(x.f[,varied.probes],dim = c(dim(x.f)[1],1)))	
	if (numPr>1)
		x.f<-x.f[,varied.probes]
	
	block<-paste("BlkCnv",1:(dim(x.f)[2]),sep="") 
	block<-paste(block,cnv.blocks$chr,sep="Chr")

	cnv.blocks<-data.frame(block=block,cnv.blocks)

	names(x.f)<-as.vector(cnv.blocks[,1])
	ans<-x.f

	attr(ans,"cnv.blocks")<-cnv.blocks
      }
	return(ans)
}
	



