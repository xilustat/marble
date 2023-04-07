
GxESelection.NonSparse=function(obj, prob=0.95){
  
  BI=obj$burn.in
  
  GS.beta = obj$posterior$GS.G
  GS.beta = GS.beta[-c(1:BI),]
  GS.eta = obj$posterior$GS.GXE
  GS.eta = GS.eta[-c(1:BI),]
  
  coef_G = c(); coef_GXE = c()
  
  for(j in 1:ncol(GS.beta))
  {
    t12 = GS.beta[,j]
    coef_G = c(coef_G,quanfun(t12,prob))
  }
  
  for(j in 1:ncol(GS.eta))
  {
    t13 = GS.eta[,j]
    coef_GXE = c(coef_GXE,quanfun(t13,prob))
  }
  
  
  
  names(coef_G) = names(obj$coefficient$G)
  names(coef_GXE) = names(obj$coefficient$GE)
  
  
  sel = list(Main.G=coef_G, GxE=coef_GXE)
  method = paste(prob*100,"% credible interval", sep = "")
  
  out = list(method=method, effects=sel)
  class(out) = "GxESelection"
  out
  
}
