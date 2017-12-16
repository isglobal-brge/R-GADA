dudi.acm.fortran<-function (df, row.w = rep(1, nrow(df)), scannf = TRUE, nf = 2) 
{
    if (!all(unlist(lapply(df, is.factor)))) 
        stop("All variables must be factors")
    X <- acm.disjonctif(df)

    lig <- nrow(X)
    col <- ncol(X)
    var <- ncol(df)
	
	
    if (length(row.w) != lig) 
        stop("Non convenient row weights")
    if (any(row.w < 0)) 
        stop("row weight < 0")
    row.w <- row.w/sum(row.w)
    col.w <- apply(X, 2, function(x) sum(x * row.w))
    if (any(col.w == 0)) 
        stop("One category with null weight")
    X <- t(t(X)/col.w) - 1
    col.w <- col.w/var
    X <- as.dudi(data.frame(X), col.w, row.w, scannf = scannf, 
        nf = nf, call = match.call(), type = "acm")
	nprobe<-ncol(df)
	nsub<-nrow(df)
	mat.f<-matrix(as.numeric(as.matrix(df)),ncol=nprobe,nrow=nsub)
	nlev<- nlev<-max(unlist(lapply(df, function(x) length(levels(x)))))
	lev<-1:nlev
	lw<-X$lw
	l1<-as.matrix(X$l1)
	rcor<-matrix(0,nrow=nprobe,ncol=nf)
	
	ans<-.Fortran("floc",
						as.double(lw),
						as.double(l1),
						as.double(mat.f),
						as.double(lev),
						as.integer(nf),
						as.integer(nprobe),
						as.integer(nsub),
						as.integer(nlev),
						rcor=as.double(rcor),
                                                PACKAGE="gada")

    rcor <- data.frame(matrix(ans$rcor,nrow=nprobe,ncol=nf))
    row.names(rcor) <- names(df)
    names(rcor) <- names(X$l1)
    X$cr <- rcor
    return(X)
}
