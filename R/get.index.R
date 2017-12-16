#get_index
#This function transforms the position of the gen porbe into an intiger index
#i.e. its position in the genetic matrix
#
#get_index -integer; index
#pos -position of the probe
#probe_pos -list of all probe positions

get.index<-function (pos,probe.pos) 
{
	index.vect<- 1:length(probe.pos)
	gen.ind<-vector(mode="numeric", length=length(pos))

	for (i in 1:length(pos))
	{
	#sometimes two SNPs have the same location! pick up the first event.
	    p<-(1:length(probe.pos))[probe.pos==pos[i]][1]
		gen.ind[i]<-index.vect[p]
	}

	ans<-sort(gen.ind)
}