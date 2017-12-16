# ----- Define a function for plotting a matrix ----- #

matrixPlot <- function(x, fac=TRUE, ...){

  #if matrix of factors convert to numeric
     if(fac)
      x<-apply(x,c(1,2),as.numeric)

     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()

  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }

# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 ColorLevels <- 0:6

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 mycol<-grey(seq(1, 0, length=7))

 image(1:length(xLabels), 1:length(yLabels), t(x), col=mycol, xlab="", ylab="Subjects", axes=FALSE,main="CNVs")

 if( !is.null(title) ){
    title(main=title)
 }

 axis(BELOW<-1, at=1:length(xLabels), labels=FALSE, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=FALSE, las= HORIZONTAL<-1, cex.axis=0.7)

 layout(1)

legend("topright", c("loss","normal", "gain"), col=mycol[c(1,3,7)], pch=15)
}
