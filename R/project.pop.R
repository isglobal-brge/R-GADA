project.pop<-function(dd,factor=FALSE,norm=FALSE, verbose=TRUE)
{
        if (factor)
         {
		   if(verbose)
		     cat("  projecting factors\n")
			 
      	   y<-as.matrix(dd$fa)
           cnt<-1:length(dd$fa[,1])
         }
        else
		 {
		   if(verbose)
		      cat("  projecting cosines\n")
			  
      	   y<-as.matrix(dd$va)
           cnt<-1:length(dd$va[,1])
         }

	#get centroids
	x<-as.matrix(dd$gc)

	ans<-list()
	
	if(norm)
	{
		L1<-diag(y %*% t(y))
		va<-as.vector(L1)

		vec.or<-order(-va)
		va.or<-va[vec.or]
		pr<-as.character(rownames(dd$va))
		pr.or<-pr[vec.or]
						

		propop<-data.frame(probe=pr,va=va)
		ans[[1]]<-propop


		for(p in 2:dim(x)[1])
		{
		ans[[p]]<-ans[[1]]
		}
		
	}
	else
	{
	#loop over number of populations
	for(p in 1:dim(x)[1])
	{
		#project
		L1<-(y %*% x[p,])/sqrt(as.vector(x[p,] %*% x[p,]))
		va<-as.vector(L1)

		vec.or<-order(-va)
		va.or<-va[vec.or]
		pr<-as.character(rownames(dd$va))
		pr.or<-pr[vec.or]
						

		propop<-data.frame(probe=pr,va=va)
		ans[[p]]<-propop
	

	}

	}
	
	nam<-rownames(x)
	names(ans)<-nam
	attr(ans,"centroid")<-dd$gc
	ans	
	
}
	
