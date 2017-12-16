collapseInfo<-function(x)
 {
  K<-x[[1]]$K
  sigma2<-x[[1]]$sigma2
  Iext<-x[[1]]$Iext
  Wext<-x[[1]]$Wext
  numit<-x[[1]]$numit
  tol<-x[[1]]$tol
  for (i in 2:length(x))
   {
    K<-c(K,x[[i]]$K)
    sigma2<- c(sigma2,x[[i]]$sigma2)
    Iext<- c(Iext,x[[i]]$Iext)
    Wext<- c(Wext,x[[i]]$Wext)
    numit<- c(numit,x[[i]]$numit)
    tol<- c(tol,x[[i]]$tol)
   }
  ans<-list(K=K,Iext=Iext,Wext=Wext,sigma2=sigma2,numit=numit,tol=tol) 
 }
