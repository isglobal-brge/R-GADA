plot.summaryParGADA<-function(x, chr, max.number.cnv=100, length.base=c(500,1e6), show.ind=FALSE, label.ind, chromosomes=c(1:22,"X","XY","Y"), ...)
 {
  
  require(plotrix)
  data(genomicInfo)
  
  if (!inherits(x, "summaryParGADA"))
   stop("object must be of class 'summaryParGADA' ")

  mypalette.blue<-c("#EFF3FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594")
  mypalette.red<-c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")

  if (chr!=0)
   load("SBL/gen.info.Rdata")

  if (chr!=0)
   selec<-c(1:length(chromosomes))[chromosomes==chr]
  else
   selec<-1

  out<-x[[selec]]
  if (chr!=0)
   {
    select<-gen.info$chr==chr
    pos<-gen.info[select,3]
    pos.all<-gen.info[,3]
   }
  else
   {
    Info<-attr(x, "Info")
    nProbes<-max(Info[[1]]$EndProbe)
    pos<-1:nProbes
    pos.all<-1:nProbes
   }

  limits<-range(pos, na.rm=TRUE, finite=TRUE)
  prob.ini<-pos[1]

  

  old.mar <- par("mar")
  old.mfrow <- par("mfrow")
  on.exit(par(mar = old.mar, mfrow = old.mfrow))

  m <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
  layout(m, heights = c(0.2, 0.8))


  nSamples<-attr(x,"Samples")[2]
  if (chr!=0)
   altered<-countAltered.i(chromosomes[selec], x, max.number.cnv, length.base, gen.info, chromosomes)
  else
   {
    altered<-countAltered.i(0, x, max.number.cnv, length.base, gen.info, 0)
    if (NROW(altered$gains)>1)
     names(altered$gains)<-c("pos","Freq")
    if (NROW(altered$losses)>1)
     names(altered$losses)<-c("pos","Freq")
   }
  
  gains<-altered$gains
  losses<-altered$losses  

#
# plot #1
#
  par(mar = c(0, old.mar[2], old.mar[3], old.mar[4]))
  plot(0,0,type="n",axes=FALSE,col="blue",ylab="% Samples",ylim=c(0,1),cex.lab=0.75,xlim=limits)

# draw percentage of losses
  if (length(losses)>0)
   {
    pos.losses<-losses$pos
    segments(pos.losses,0,pos.losses,losses$Freq/nSamples,col="blue")
   }

# draw percentage of gains
  if (length(gains)>0)
   {
    pos.gains<-gains$pos
    segments(pos.gains,0,pos.gains,gains$Freq/nSamples,col="red")
   }


  cc<-c(0,.25,.5,.75,1)
  axis(2,at=cc,label=c(".0",".25",".50",".75","1"),cex.axis=0.8,las=1)
  title(paste("Chromosome",chr),line=1.5)
  title(paste("(",nSamples," samples)",sep=""),line=0.5,cex.main=0.8)

  segments(min(pos, na.rm=TRUE), 0, max(pos, na.rm=TRUE), 0)
  for (i in 2:5)
   {
    segments(min(pos, na.rm=TRUE),cc[i],max(pos, na.rm=TRUE),cc[i],col="gray75")
   }

#
# ----------------------------------------------------------------
#
  draw.segments<-function(i,x,thresholds,palette,nSamples,limits,show.ind,label.ind,max.number.cnv,length.base)
   {
     seg<-x[[i]]
     mycol<-as.character(as.factor(cut(seg[,4],thresholds,labels=palette)))
     J<-nrow(seg)
     control<-ifelse(nrow(seg)<max.number.cnv, TRUE, FALSE)
     if (J>0 & control)
      {
       for (j in 1:J)
        {
          if (i > 1 & i < nSamples)
           {
            ini <- i - .5
            end <- i + .5
            if (show.ind)
             {
              segments(limits[1],ini,limits[2],ini,lty=2,col="gray50") 
              segments(limits[1],end,limits[2],end,lty=2,col="gray50")
             }
          }
          else if (i==1)
           {
            ini <- i 
            end <- i + .5
           }
          else
           {
            ini <- i - .5
            end <- i
           }

          ll<-seg[j,2]-seg[j,1] 
          if (length.base[1]<ll & ll<length.base[2]) 
           {
             polygon(c(seg[j,1],seg[j,2],seg[j,2],seg[j,1]),c(ini,ini,end,end),
                col=mycol[j],border=mycol[j])
           } 

          text(limits[1],i,label.ind[i],cex=0.5,adj=1)


        }
     }

   }
  
#
# plot #2
#
  par(mar = c(old.mar[1], old.mar[2], 0, old.mar[4]))
  size.ok<-nSamples/50
  drawChromosome(chr, size=size.ok, ylim=c(1,nSamples), limits=limits, ylab="Individuals")


  thresholds<-c(-Inf,0,Inf)
  palette<-c(mypalette.blue[5],mypalette.red[5])

  if (missing(label.ind) & !show.ind)
   label.ind<-rep("",nSamples)
  
  if (missing(label.ind) & show.ind)
   label.ind<-labels(x)


  invisible(lapply(1:nSamples, draw.segments, x=out, thresholds=thresholds, palette=palette, nSamples=nSamples, limits=limits, show.ind=show.ind, label.ind=label.ind, max.number.cnv=max.number.cnv, length.base=length.base) )
 
 }

