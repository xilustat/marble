
BL_SS <- function(X, Y, E, clin, max.steps, debugging=FALSE)
{
  dat = DataMatrix(X, Y, E, clin, intercept=TRUE, debugging=FALSE)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y;
  n = dat$n; p = dat$p; q=ncol(c)
  env = dat$env
  
  G.names = dat$G.names
  E.names = dat$E.names
  clin.names = dat$clin.names
  GXE.names = dat$GXE.names
  
  #max.steps=10000
  hatAlpha = rep(1,env); hatb = rep(1,q); hatEta= rep(1,env); hatBeta=1
  invSigAlpha0= diag(rep(1,env)); invSigb0 = diag(rep(1,q))
  hatInvTauSq1= 1; hatInvTauSq2 = rep(1,env)
  sg1=1; sg2 = rep(1,env); hatLambdaSqStar1=1; hatLambdaSqStar2=1; hatSigmaSq=1
  aStar=1; bStar=1; alpha=1; gamma=1
  progress = 0; hatPiEta=1/2; hatPiBeta=1/2
  mu0=1; nu0=1
  
  coef_G_BLSS = c(); coef_GXE_BLSS = c(); 
  coef_G_SS = c(); coef_GXE_SS = c(); 
  
  for(j in 1:p)
  {
    x=g[,j]
    w=xx[,((env*(j-1)+1):(j*env))]
    betas = BLSS(x,y,e,c,w,max.steps,n,hatBeta,hatEta,hatAlpha,hatb,hatInvTauSq1,hatInvTauSq2,sg1,sg2,hatPiEta,hatPiBeta,invSigAlpha0, invSigb0, hatLambdaSqStar1, hatLambdaSqStar2,hatSigmaSq, aStar, bStar, alpha, gamma,mu0,nu0, progress)
    t1=as.matrix(betas$GS.beta)
    coef_G_BLSS=cbind(coef_G_BLSS, t1)
    t12= betas$GS.SS1
    t12 = t12[seq(max.steps/2, max.steps,1)]
    q_beta = mean(t12)
    
    coef_G_SS = c(coef_G_SS,q_beta)
    
    coef_eta=betas$GS.SS2
    
    for(d in 1:env)
    {
      t2 = as.matrix(betas$GS.eta[,d])
      coef_GXE_BLSS = cbind(coef_GXE_BLSS, t2)
      t13=as.matrix(coef_eta[,d])
      t13 = t13[seq(max.steps/2, max.steps,1),]
      q_eta = mean(t13)
      
      coef_GXE_SS = c(coef_GXE_SS,q_eta)
    }
  
  }
  
  BI=max.steps/2
  coeff.E = apply(betas$GS.alpha[-(1:BI),,drop=FALSE], 2, stats::median);
  names(coeff.E) = E.names;
  coeff.clin = apply(betas$GS.b[-(1:BI),,drop=FALSE], 2, stats::median); 
  names(coeff.clin) = c(1,clin.names);
  
  coeff.G = apply(coef_G_BLSS[-(1:BI),,drop=FALSE], 2, stats::median); 
  names(coeff.G) = G.names;
  coeff.GE = apply(coef_GXE_BLSS[-(1:BI),], 2, stats::median); 
  names(coeff.GE) = GXE.names;
  
  coef_hat = c(coef_G_SS, coef_GXE_SS);names(coef_hat)=c(G.names,GXE.names)
  rank.GXE = sort(coef_hat,decreasing = TRUE)
  
  coefficient = list(clin=coeff.clin, E=coeff.E, GE=coeff.GE, G=coeff.G)
  out = list(GS.E=betas$GS.alpha, GS.C=betas$GS.b, GS.G=coef_G_BLSS, GS.GXE=coef_GXE_BLSS)
  fit = list(posterior = out, coefficient=coefficient, ranklist=rank.GXE, burn.in = BI, iterations=max.steps, design=list(g,xx=xx, CLC=c, E=e))
  
  fit
  
}