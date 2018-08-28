summary.BackwardElimination <-function(object, which="both", T, BaseAmp, print=TRUE, ...)
 {
  which.ok<-match(which,c("Gains","Loses","both"),nomatch=0)
  if (which.ok==0)
   stop ("which should be 'Gains', 'Loses' or 'both' ")
  
  gen.info <- attr(object,"gen.info")

#RPR Pulling out the segments, and segments amplitudes
  if (!is.null(gen.info))
   {
    k <- attr(object, "chr")
    Segments<- WextIextToSegments(object[[1]])
    Segments$chromosome <- attr(object,"chr")[[1]]
    for (i in 2:length(k))
     {
      temp<-WextIextToSegments(object[[i]]) 
      temp$chromosome <- attr(object,"chr")[[i]]
      temp[,1]<-temp[,1]+Segments[nrow(Segments),2]  
      temp[,2]<-temp[,2]+Segments[nrow(Segments),2]
      Segments<-rbind(Segments,temp)
     }
    }
  else 
   {
    Segments<-WextIextToSegments(object)
    Segments$chromosome <- 0 #attr(object,"chr")[[i]] #FIXME If no genome.info.... I don't know the chr... maybe we can set 0 as a mark... 
   }

   K <- nrow(Segments)
 
    sigma2<-object$sigma2
# Estimation of the reference level as the median across all probes
   if (missing(BaseAmp)) 
   {


# JRG Oct'09
    if (!is.null(gen.info))
     {
       Segments.autosomes<-Segments[Segments$chromosome%in%c(1:22),]
       Segments.X<-Segments[Segments$chromosome%in%c("X"),]
       K.X<-nrow(Segments.X)
       Segments.Y<-Segments[Segments$chromosome%in%c("Y"),]
       K.Y<-nrow(Segments.Y)
     }

    else
     {
       Segments.autosomes<-Segments
     }
  

     K.autosomes<-nrow(Segments.autosomes)


# For autosomes 
     aux<-.C("RcallCompAmpMedianMethod",
            SegLen=as.integer(Segments.autosomes$LenProbe),
            SegAmp=as.double(Segments.autosomes$MeanAmp),
            K=as.integer(K.autosomes),
            BaseAmp=double(1),PACKAGE="gada")

     BaseAmp <- ifelse(is.na(aux$BaseAmp),0,aux$BaseAmp)


# For X and Y chromosomes
    if (!is.null(gen.info))
     {
       aux.X<-.C("RcallCompAmpMedianMethod",
              SegLen=as.integer(Segments.X$LenProbe),
              SegAmp=as.double(Segments.X$MeanAmp),
              K=as.integer(K.X),
              BaseAmp=double(1),PACKAGE="gada")

       BaseAmp.X <- ifelse(is.na(aux.X$BaseAmp),0,aux.X$BaseAmp)


       aux.Y<-.C("RcallCompAmpMedianMethod",
              SegLen=as.integer(Segments.Y$LenProbe),
              SegAmp=as.double(Segments.Y$MeanAmp),
              K=as.integer(K.Y),
              BaseAmp=double(1),PACKAGE="gada")

       BaseAmp.Y <- ifelse(is.na(aux.Y$BaseAmp),0,aux.Y$BaseAmp)

       BaseAmp.all<-c(BaseAmp,BaseAmp.X,BaseAmp.Y)
      }
     
     else
      {
       BaseAmp.all<-c(BaseAmp, NA, NA)
      }

   }

   else
    {
      if (length(BaseAmp)!=3)
        stop("'BaseAmp' should have three components: autosomes, X, Y")

      BaseAmp.all<-BaseAmp
    }



  if (!is.null(gen.info))
  {
# if T missing use x$T
   if (missing(T))
     T <- object[[1]]$T
   sigma2 <- object[[1]]$sigma2
   MinSegLen <- object[[1]]$MinSegLen
  }

  else
  {
# if T missing use x$T
   if (missing(T))
     T <- object$T
   sigma2 <- object$sigma2
   MinSegLen <- object$MinSegLen
  }


   aux<-.C("RcallClassifySegments",
            SegLen=as.integer(Segments.autosomes$LenProbe),            
            SegAmp=as.double(Segments.autosomes$MeanAmp),
            SegState=double(K.autosomes),
            K=as.integer(K.autosomes),
            BaseAmp=as.double(BaseAmp.all[1]),
            sigma2=as.double(sigma2),
            T=as.double(T),PACKAGE="gada")                   

   Segments.autosomes$State <- aux$SegState




  if (!is.null(gen.info)) # JRG Oct'09
    {

      aux.X<-.C("RcallClassifySegments",  
               SegLen=as.integer(Segments.X$LenProbe),            
               SegAmp=as.double(Segments.X$MeanAmp),
               SegState=double(K.X),
               K=as.integer(K.X),
               BaseAmp=as.double(BaseAmp.all[2]),
               sigma2=as.double(sigma2),
               T=as.double(T),PACKAGE="gada")                   


      aux.Y<-.C("RcallClassifySegments",
               SegLen=as.integer(Segments.Y$LenProbe),            
               SegAmp=as.double(Segments.Y$MeanAmp),
               SegState=double(K.Y),
               K=as.integer(K.Y),
               BaseAmp=as.double(BaseAmp.all[3]),
               sigma2=as.double(sigma2),
               T=as.double(T),PACKAGE="gada")                   
  
      Segments$State <- c(aux$SegState, aux.X$SegState, aux.Y$SegState)
    }

  else
   {
     Segments$State <- aux$SegState
   }

   
#
# Get genomic information
#
  if (!is.null(gen.info)) # JRG Oct'09
    {
      ini <- as.character(Segments[,1]) # changed 08/2018 bad annot
      end <- as.character(Segments[,2]) # changed 08/2018 bad annot
      Segments[,1] <- gen.info[ini,3]
      Segments[,2] <- gen.info[end,3]
    }
  else
   {
     ini<-Segments[,1]
     end<-Segments[,2]
   }


#
# Normalize Segments JRG
#

# Autosomes
  if (!is.null(gen.info))
    {
      o<-Segments[,5]%in%c(1:22)
      Segments[o,4]<-Segments[o,4]-BaseAmp.all[1]

if (F) # BaseAmp.all is not used to normalize, only to classify segments.
 {
  # Chr X,Y
   o<-Segments[,5]%in%c("X")
   Segments[o,4]<-Segments[o,4]-BaseAmp.all[2]
   o<-Segments[,5]%in%c("Y")
   Segments[o,4]<-Segments[o,4]-BaseAmp.all[3]
 } 

    }

  else

    {
      Segments[,4]<-Segments[,4]-BaseAmp.all[1]

    }

  

   
  if (print)
   {
    cat("---------------------------------------- \n") 
    cat("Sparse Bayesian Learnig (SBL) algorithm \n")
    cat("Backward Elimination procedure with T=",T," and minimun length size=",MinSegLen, "\n",sep="")
    cat(" Number of segments = ",K,"\n")
    cat(" Base Amplitude of copy number 2: chr 1:22:", round(BaseAmp.all[1],4), ", X=", round(BaseAmp.all[2],4), ", Y=",round(BaseAmp.all[3],4),"\n", sep="")
    
    if (which.ok==1)
     {
      cat("Gains with respect Base Amplitude \n")
      cat("---------------------------------------- \n") 
      print(Segments[Segments$State==1,])
     }

    if (which.ok==2)
     {
      cat("Loss with respect Base Amplitude \n")
      cat("---------------------------------------- \n") 
      print(Segments[Segments$State==-1,])
     }

    if (which.ok==3)
     {
      cat(" Gains (1) and Loses (-1) with respect Base Amplitude \n")
      cat("---------------------------------------- \n") 
      print(Segments[Segments$State%in%c(1,-1),])
     }
   }
  else
   {
    cat("---------------------------------------- \n") 
    cat("Sparse Bayesian Learnig (SBL) algorithm \n")
    cat("Backward Elimination procedure with T=",T," and minimun length size=",MinSegLen, "\n",sep="")
    cat(" Number of segments = ",K,"\n")
    cat(" Base Amplitude of copy number 2: chr 1:22:", round(BaseAmp.all[1],4), ", X=", round(BaseAmp.all[2],4), ", Y=",round(BaseAmp.all[3],4),"\n", sep="")
   }
  
 
 class(Segments)<-c("data.frame","summary.BackwardElimination")
 names(BaseAmp.all)<-c("Autosomes","X","Y")
 attr(Segments,"BaseAmp")<-BaseAmp.all
 attr(Segments,"index")<-cbind(ini,end)
 if (is.null(gen.info))
  sigma2<-object$sigma2
 else
  sigma2<-object[[1]]$sigma2
 attr(Segments,"sigma2")<-sigma2
 invisible(Segments)

 }



