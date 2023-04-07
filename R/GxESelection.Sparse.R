
GxESelection.Sparse=function(obj){
  BI=obj$burn.in
  
  GS.beta = obj$posterior$GS.G
  GS.beta = GS.beta[-c(1:BI),]
  GS.eta = obj$posterior$GS.GXE
  GS.eta = GS.eta[-c(1:BI),]
  
  coef_G = c(); coef_GXE = c()
  
  for(j in 1:ncol(GS.beta))
  {
    t12 = GS.beta[,j]
    t12[t12!=0]=1; t12[t12==0]=0
    q_beta = mpm(t12)
    coef_G = c(coef_G,q_beta)
  }
  
  for(j in 1:ncol(GS.eta))
  {
    t13 = GS.eta[,j]
    t13[t13!=0]=1; t13[t13==0]=0
    q_eta = mpm(t13)
    coef_GXE = c(coef_GXE,q_eta)
  }
  
  
  
  names(coef_G) = names(obj$coefficient$G)
  names(coef_GXE) = names(obj$coefficient$GE)
  
  
  sel = list(Main.G=coef_G, GxE=coef_GXE)
  method = paste("posterior inclusion Proportion", sep = "")
  
  out = list(method=method, effects=sel)
  class(out) = "GxESelection"
  out
}

