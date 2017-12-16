getChromosomeDat<-function(object, ch, gg, verbose=TRUE)
{
	cat("Getting data from chromosome:", ch, "\n")

	probe.pos<-gg[gg$chr==ch,]$pos
	names(probe.pos)<-gg[gg$chr==ch,]$probe

	cnv.data<-lapply(1:length(object), function(x) object[[x]][[ch]])
	cnv.data<-unlist(cnv.data,recursive=FALSE)
	
	mat<-get.matrix(cnv.data, probe.pos, verbose)
}

