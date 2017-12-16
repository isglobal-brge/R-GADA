
plotWG<-function(x, min.percentage=0.05, max.number.cnv=100, length.base, chr=c(1:22,"X","Y"))
 {

   if (!inherits(x, "summaryParGADA")) 
    stop("object must be of class 'summaryParGADA'")

 
   if (missing(length.base))
    length.base<-attr(x,"length.base")

   require(plotrix)   ##RPR interferes with Aroma.Affymetrix...
                      ##JRG this is why it is here
   data(genomicInfo)

   if (!"gen.info"%in%ls(pos=.GlobalEnv))
    {
     Info<-attr(x,"Info")
     load(paste(Info,"/SBL/gen.info.Rdata",sep=""))
    }

   col.legend<-c("red","blue") 
   chrs <- gen.info[, 2]
   limits <- c(min(gen.info[, 3]), max(gen.info[, 3]))
  

   n.chr <- length(chr)

   old.mfrow <- par("mfrow")
   old.mar <- par("mar")
   old.xpd <- par("xpd")
   on.exit(par(mfrow = old.mfrow, mar = old.mar))
   par(mfrow = c(n.chr+2, 1))
   par(mar = c(0.5, 5, 0.2, 3))
   par(xpd = TRUE)
   plot(c(1:3), rep(1, 2, 3), axes = FALSE, xlab = "", ylab = "", type = "n")
   legend(1, 1, c("Gains", "Losses"), col = col.legend, 
            pt.bg = col.legend, pch = rep(22, 4), horiz = TRUE, 
            cex = 1, bty = "n", pt.cex = 1.6, yjust = 0.5)

   nSamples<-attr(x,"Samples")[2]

   text(2,1,paste("CNV frequency summary (",nSamples," samples)",sep=""),font=2)
         
  for (i in 1:n.chr) 
   {
      select <- gen.info[, 2] == chr[i]
      pos.chr <- gen.info[select, 3]  


      drawChromosome(i,0.5,print.names=FALSE,ylim=c(-0.5,1),limits=limits,ylab="")

            ss<-c(0,0.25,0.5,0.75,1)  
            for (j in 1:5)
             {  
               segments(limits[1], ss[j], max(pos.chr), ss[j], col="gray60")
             }



# count gains and losses
            altered<-countAltered.i(chr[i], x, max.number.cnv, length.base, gen.info)

            if (!is.null(altered$gains))
             {
# draw gains
              gains<-altered$gains
              pos.gains<-gains$pos
              hh2<-gains$Freq/nSamples
              hh2[hh2<min.percentage]<-NA
              if (length(pos.gains)>0)
               segments(pos.gains, 0, pos.gains, hh2, col=col.legend[1])
             }


            if (!is.null(altered$losses))
             {
# draw losses
              losses<-altered$losses
              pos.losses<-losses$pos
              hh<-losses$Freq/nSamples
              hh[hh<min.percentage]<-NA
              if (length(pos.losses)>0)
               segments(pos.losses, 0, pos.losses, hh, col=col.legend[2])
             }


            text(par("usr")[1], 0, chr[i], cex = 1, adj = 0.5)


            text(par("usr")[1], 0, chr[i], cex = 1, adj = 0.5)
        }
 
    
        plot(limits, c(1, 1), axes = FALSE, xlab = "", ylab = "", 
            type = "n", xlim = limits)
        text(limits[1], 1, limits[1], cex = 1, adj = 0)
        text(limits[2], 1, limits[2], cex = 1, adj = 1)
        text((limits[2] - limits[1])/2, 1, "Genomic Position", 
            cex = 1, adj = 0, font = 2)
      
    }



