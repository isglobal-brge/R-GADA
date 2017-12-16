discriminCNV<-function(mat.f, pop.cla, select=NULL, type="acm", verbose=TRUE, ...)
{

   proj<-getProj(mat.f, pop.cla, select=select, saveInfo=TRUE, type=type,verbose,...)
   cen<-attr(proj, "centroid")

   if(!all(unlist(lapply(proj,function(x) names(x)==c("probe", "va")))))
		stop("projection over populations expected")
   if (!is.data.frame(mat.f)) 
    stop("data.frame expected")

   if (!is.factor(pop.cla))
    stop("pop must be a factor")
		
   ans<-list(proj=proj, pop=pop.cla, cen=cen, mat.f=mat.f)	
	
   class(ans)<-"cnv.discrimin"
	
   return(ans)
}
