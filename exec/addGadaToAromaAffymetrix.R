####################################################################################
##  Adding GADA in Aroma.Affymetrix
####################################################################################

## 1) Create a subclass GadaModel of the abstract class
## CopyNumberSegmentationModel, to set up the model.

setConstructorS3("GadaModel", function(cesTuple=NULL, ...,
  aAlpha=0.2,estim.sigma2=TRUE,T=3,MinSegLen=3) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cesTuple)) {
    require('gada') || throw("Package not loaded: GADA");
  }

  # We pick up the parameters for the GADA model
  GadaArgs=list(aAlpha=aAlpha,estim.sigma2=estim.sigma2,T=T,MinSegLen=MinSegLen)
  
  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "GadaModel",
         GadaArgs=GadaArgs)
})

## 2) provide a method getAsteriskTags() that returns the
## tag(s) replacing any tag "*":
setMethodS3("getAsteriskTags", "GadaModel", function(this, collapse=NULL, ...) {
  tags <- "GADA";

  # Add class-specific tags
  GadaArgs <- this$GadaArgs;
  tags <- c(tags, paste("a",GadaArgs$aAlpha,sep=""));
  
  if (isPaired(this))
    tags <- c(tags, "paired");

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)

## Now we can already try assuming we have a
## CnChipEffectSet object named 'ces' :
## > gada <- GadaModel(ces);
## > print(gada);
## > print(getPath(gada)); 

## Note how the "root" path of your GadaModel starts with 'gadaData/'.
## That is the output path where the data is stored containing the model
## estimates/results.  If we download the source of aroma.affymetrix 
## we can see how this is done for CbsModel in source file R/CbsModel.R.


#####################################################################################

## 3) fitOne() method to fit the model

## We can now even try the following:
## > fit(gada, arrays=1:2, chromosomes=22, verbose=-20);
## We can see it that the framework will load the chip effects for the
## chromosomes of interest, calculate the raw CN (log2 ratios) and try to
## call fitOne() on our GadaModel object.  We get an
## error because we still haven't defined fitOne() for the GadaModel
## class.  This is the method that will take the data and do the fitting
## etc and return the results. 

## fitOne() should take an argument 'data' which is a Jx6 matrix where
## the number of rows J is the number of loci, and the six columns are x,
## M, sdTheta, sdM, chipType, and unit.  To understand what this object
## is, try:
## > files <- getMatrixChipEffectFiles(gada, array=1);
## > ceList <- files[,"test"];
## > refList <- files[,"reference"];
## > data <- getRawCnData(gada, ceList=ceList, refList=refList, chromosome=22, verbose=-10);
## > str(data)
##  num [1:74, 1:6] 15685581 16604328 16958128 17506879 18008688 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : NULL
##   ..$ : chr [1:6] "x" "M" "sdTheta" "sdM" ...
## > plot(data[,c("x","M")], ylim=c(-3,3))

##   'x' is the genomic position,
##   'M' is the log2 ratio
##   'sdM' ('sdTheta') is the estimated standard of each M (theta), and
##   'chiptype' and 'unit' is the chip type and the CDF unit index for the particular locus.
##       'chiptype' is useful if data from two or more chip types are to be segmented.
## Assume you have an R function fitGada(x,M) that takes positions 'x'
## and log-ratios 'M' and identifies CN regions returned in some
## structure.  For simplicity, assume the return structure is of class
## "GadaFit".  Then fitOne() should look something like:
##   setMethodS3("fitOne", "GadaModel", function(this, data, ..., verbose=FALSE) {
##       fit <- fitGada(x=data[,"x"], M=data[,"M"]);
##         # In case your fitGada() does not return an object with class, then do
##         # class(fit) <- "GadaFit"
##         # In case your fitGada() does not return the raw CNs, add them here
##         # fit$data <- data;
##         fit;
##     }, private=TRUE) # fitOne()

## You probably want to add verbose output etc.  See how this is done for
## See the GladModel or CbsModel (see R/CbsModel.fitOne.R) for another example.

setMethodS3("fitOne", "GadaModel", function(this, data, chromosome, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Fitting GADA");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract arguments for GADA functions SBL() and BackwardElimination()
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  args <- this$GadaArgs; ## If passed through GadaModel
  #args <- list(...);    ## If passed through fit argument list
  keep <- (names(args) %in% names(formals(gada::SBL))); 
  fitSblArgs <- args[keep];
  keep <- (names(args) %in% names(formals(gada::BackwardElimination))); 
  fitBackElimArgs <- args[keep];

  
  verbose && enter(verbose, "Setting up GADA data structure");
  nbrOfUnits <- nrow(data);
  chipTypes <- getChipTypes(this);
  ## Setting up the GADAobject mydata  manually instead of setupGADA...
  mydata<- list();
  mydata$log.ratio<-data[,"M"];
  verbose && cat(verbose, "The  ",sum(is.na(mydata$log.ratio))," NAs on the data are interpreted as log2.ratio=0"); #Move down
  attr(mydata,"gen.info") <- NULL; 
  mygenome.info <- data.frame(probe=data[,"unit"], #Consider add chipType?
                              chr=rep(chromosome,nbrOfUnits),
                              pos=data[,"x"],
                              stringAsFactors=FALSE);                                             
  attr(mydata, "type") <-"log.ratio";
  class(mydata) <- "setupGADA";
  ## Done with the object
  
  verbose && str(verbose, mydata);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling SBL()"); 
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Chip types: ", paste(chipTypes, collapse=", "));
  verbose && cat(verbose, "Total number of units: ", nbrOfUnits); 

  rm(data);
# Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

# SBL() writes to stdout only with verbose; capture it and send it to the verbose object.  
   stdout <- capture.output({ 
        SBLfit <- do.call("SBL",c(list(mydata),fitSblArgs));  
   })
   stdout <- paste(stdout, collapse="\n");
   verbose && cat(verbose, stdout);
        
   verbose && cat(verbose,print(SBLfit)); 
   verbose && exit(verbose);
  
  verbose && enter(verbose, "Calling BackwardElimination()");
  stdout <- capture.output({ 
         fit <- do.call("BackwardElimination",c(list(SBLfit),fitBackElimArgs));  
   })
   stdout <- paste(stdout, collapse="\n");
   verbose && cat(verbose, stdout); 
  
  verbose && cat(verbose,print(fit)); # This gives verbose info.
##  verbose && cat(verbose,summary(fit)); # This is too much detail, i.e the segments
  verbose && exit(verbose);

  verbose && exit(verbose);

  attr(fit,"mygenome.info") <- mygenome.info; #supplying probe info....
  attr(fit,"GadaArgs") <- this$GadaArgs; #supply model parameters.
  class(fit) <- "GadaFit"  #
  
  fit;  
}, private=TRUE) # fitOne()


################################################################################

## 3) Define extractRawCopyNumbers() for "GadaFit" that returns a RawCopyNumber object:
##    In theory this could be provided without the GadaFit object but this is how
##    aroma.affymetrix is designed.

setMethodS3("extractRawCopyNumbers", "GadaFit", function(object, ...) {
  genome.info <- attr(object,"mygenome.info") 
  chromosome <- unique(genome.info$chr);  # Can I use first element instead of unique?  
  RawCopyNumbers(cn=attr(object,"data"), x=genome.info$pos, chromosome=chromosome);
})

################################################################################

## 4) extractCopyNumberRegions() for "GadaFit" that returns a
## CopyNumberRegions object, e.g.

setMethodS3("extractCopyNumberRegions", "GadaFit", function(object, ...) {

  mygenome.info<-attr(object,"mygenome.info");
  GadaArgs<-attr(object,"GadaArgs");

## We can use the arguments to obtain a higher T without having to fit the model again.
## but this is only possible for a T, and MinSegLen larger than than that on the stored model  
  args <- list(...);
##  keep <- (names(args) %in% names(formals(gada::SBL))); # We cannot change this parameters of the fit
##  fitSblArgs <- args[keep];
##  if length(fitSblArgs)>0 we need to repeat the fit !!  

## We can change the Backward Elimination stoping criteria to a higher value without repeating the fit.
  keep <- (names(args) %in% names(formals(gada::BackwardElimination))); 
  fitBackElimArgs <- args[keep];

  if (length(fitBackElimArgs)>0){ # Call BackwardElimination on object again
    rep=1;
    if (!("T" %in% names(fitBackElimArgs)))
        fitBackElimArgs$T=GadaArgs$T;
    if (!("MinSegLen" %in% names(fitBackElimArgs)))
        fitBackElimArgs$MinSegLen=GadaArgs$MinSegLen;
    if (fitBackElimArgs$T<GadaArgs$T){
      throw("Current Fit uses a large T and cannot be adjsted to a lower value, you need to repeat fit with a lower value");
      rep=0;
    }
    if (fitBackElimArgs$MinSegLen<GadaArgs$MinSegLen){
      throw("Current Fit uses a large MinSegLen and cannot be adjusted to a lower value, you need to repeat fit with a lower value");
      rep=0;
    }
    if (rep==1){
      verbose && enter(verbose, "Repeating BackwardElimination()");
      verbose && cat(verbose,str(fitBackElimArgs));
      class(object)<-"SBL";
      stdout <- capture.output({ 
        object <- do.call("BackwardElimination",c(list(object),fitBackElimArgs));  
      })
      class(object)<-"GadaFit";
      stdout <- paste(stdout, collapse="\n");
      verbose && cat(verbose, stdout);
    }  
  }

  Segments<-WextIextToSegments(object);  
        
  CopyNumberRegions(
     chromosome= mygenome.info$chr[1], 
     start=mygenome.info$pos[Segments$IniProbe], 
     stop=mygenome.info$pos[Segments$EndProbe], 
     mean=Segments$MeanAmp,
     count=Segments$LenProbe  ##Remember I can add $call at the end! using a simple threshold
  );

})

