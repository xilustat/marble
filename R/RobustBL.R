
RobustBL <- function(X, Y, E, clin, max.steps, debugging=FALSE)
{
  dat = DataMatrix(X, Y, E, clin, intercept=TRUE, debugging=FALSE)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y;
  n = dat$n; p = dat$p; q=ncol(c)
  env = dat$env
  
  G.names = dat$G.names
  E.names = dat$E.names
  clin.names = dat$clin.names
  GXE.names = dat$GXE.names
  
  hatAlpha = rep(1,env); hatb = rep(1,q); hatEta= rep(1,env); hatBeta=1; hatTau=1; hatV = rep(1,n)
  invSigAlpha0= diag(rep(1,env)); invSigb0 = diag(rep(1,q))
  hatSg1=1; hatSg2 = rep(1,env); hatEtaSq1=1; hatEtaSq2=1
  r1=1; r=1; a=1; b=1; theta=0.5
  progress = 0
  
  coef_G_RBL = c(); coef_GXE_RBL = c(); 
  
  for(j in 1:p)
  {
    
    x=g[,j]
    w=xx[,((env*(j-1)+1):(j*env))]
    betas = RBL(x,y,w,c,e,max.steps,n,hatAlpha,hatb,hatBeta,hatEta,hatTau,hatV,hatSg1,hatSg2,invSigAlpha0, invSigb0, hatEtaSq1, hatEtaSq2, theta, r1, r,a ,b, progress)
    t1=as.matrix(betas$GS.beta)
    coef_G_RBL=cbind(coef_G_RBL, t1)
    
    for(d in 1:env)
    {
      t2 = as.matrix(betas$GS.eta[,d])
      coef_GXE_RBL = cbind(coef_GXE_RBL, t2)
    }
    
    
  }
  
  
  prob_seq=seq(0.2,1,by=0.05)
  coef_G_SS_temp=c()
  
  for (k in 1:length(prob_seq)) 
  {
    beta_hat = c()
    prob=prob_seq[k]
    for (j in 1:ncol(coef_G_RBL)) 
    {
      t11 = coef_G_RBL[,j]
      beta_hat = c(beta_hat,quanfun(t1,prob))
    }
    coef_G_SS_temp=rbind(coef_G_SS_temp,beta_hat)
  }
  coef_G_SS=colMeans(coef_G_SS_temp)
  
  coef_GXE_SS_temp=c()
  for (k in 1:length(prob_seq)) 
  {
    eta_hat = c()
    prob=prob_seq[k]
    for (j in 1:ncol(coef_GXE_RBL)) 
    {
      t12 = coef_GXE_RBL[,j]
      eta_hat = c(eta_hat,quanfun(t12,prob))
    }
    coef_GXE_SS_temp=rbind(coef_GXE_SS_temp,eta_hat)
  }
  coef_GXE_SS=colMeans(coef_GXE_SS_temp)
  
  
  BI=max.steps/2
  coeff.E = apply(betas$GS.alpha[-(1:BI),,drop=FALSE], 2, stats::median);
  names(coeff.E) = E.names;
  coeff.clin = apply(betas$GS.b[-(1:BI),,drop=FALSE], 2, stats::median); 
  names(coeff.clin) = c(1,clin.names);
  
  coeff.G = apply(coef_G_RBL[-(1:BI),,drop=FALSE], 2, stats::median); 
  names(coeff.G) = G.names;
  coeff.GE = apply(coef_GXE_RBL[-(1:BI),], 2, stats::median); 
  names(coeff.GE) = GXE.names;
  
  coef_hat = c(coef_G_SS, coef_GXE_SS);names(coef_hat)=c(G.names,GXE.names)
  rank.GXE = sort(coef_hat,decreasing = TRUE)
  
  coefficient = list(clin=coeff.clin, E=coeff.E, GE=coeff.GE, G=coeff.G)
  out = list(GS.E=betas$GS.alpha, GS.C=betas$GS.b, GS.G=coef_G_RBL, GS.GXE=coef_GXE_RBL)
  fit = list(posterior = out, coefficient=coefficient,ranklist=rank.GXE, burn.in = BI, iterations=max.steps, design=list(g,xx=xx, CLC=c, E=e))
  
  fit
}
