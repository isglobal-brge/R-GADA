assocCNV.i<-function(x, formula, family=binomial, ...)
 {
   assign("x", x, envir = .GlobalEnv)
   mod0<-glm(formula, family=family, ...)
   mod1<-update(mod0, . ~ . + x)
   res<-anova(mod0, mod1, test="Chi")$P[2]
   res
 } 

