print.cnv.discrimin<-function(x, ...)
{
	cat("CNV projection \n")
	cat("class:")
    cat(class(x))
	cat("\n")
	
	cat("\n")
	sumry <- array("", c(1, 4), list(1:1, c("list", "length", "attr","content")))
    sumry[1, ] <- c("$proj", length(x$proj),"dudi.dis" ,"projection on centroids")
	rownames(sumry)<-c("")
	class(sumry) <- "table"
    print(sumry)
	cat("\n")
	
	cat("	attr: dudi.dis is class dudi")
	cat("\n")
	
	
	cat("\n")
	
	sumry <- array("", c(2, 4), list(1:2, c("data.frame", "nrow", 
        "ncol", "content")))
	sumry[1, ] <- c("$cen", nrow(x$cen), ncol(x$cen), "centroid coordinates")
	
	if(sum(dim(x$mat.f))==0)
	{
		sumry[2, ] <- c("", "","", "")
		rownames(sumry)<-c("","")
	}
	else
		sumry[2, ] <- c("$mat.f", nrow(x$mat.f), ncol(x$mat.f), "reduced matrix")
		
	class(sumry) <- "table"
    print(sumry)
    cat("\n")

	sumry <- array("", c(1, 4), list(1:1, c("factor", "length", "levels","content")))
    sumry[1, ] <- c("$pop", length(x$pop), length(levels(x$pop)), "population membership")
	rownames(sumry)<-c("")
	class(sumry) <- "table"
    print(sumry)
    cat("\n")

	
}
