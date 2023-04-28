
RobustBLSS <- function(X, Y, E, clin, max.steps, debugging=FALSE)
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
  theta=0.5; r1=1; r=1; a=1; b=1; sg1=1; sg2 = rep(1,env)
  hatPiBeta=1/2; hatPiEta=1/2
  sh0=1; sh1=1
  progress = 0
  
  coef_G_RBLSS = c(); coef_GXE_RBLSS = c(); 
  coef_G_SS = c(); coef_GXE_SS = c(); 
  
  for(j in 1:p)
  {
    
    x=g[,j]
    w=xx[,((env*(j-1)+1):(j*env))]
    
    betas = RBLSS(x,y,w,c,e,max.steps,n,hatAlpha,hatb,hatBeta,hatEta,hatTau,hatV,hatSg1,hatSg2,sg1,sg2,invSigAlpha0, invSigb0, hatEtaSq1, hatEtaSq2, theta, r1, r,a ,b ,hatPiBeta,hatPiEta,sh0,sh1, progress)
    t1=as.matrix(betas$GS.beta)
    coef_G_RBLSS=cbind(coef_G_RBLSS, t1)
    t12 = t1[seq(max.steps/2, max.steps,1)]
    t12[t12!=0]=1; t12[t12==0]=0
    q_beta = mean(t12)
    coef_G_SS = c(coef_G_SS,q_beta)
    
    coef_eta=betas$GS.eta
    for(d in 1:env)
    {
      t2 = as.matrix(betas$GS.eta[,d])
      coef_GXE_RBLSS = cbind(coef_GXE_RBLSS, t2)
      
      t13 = t2[seq(max.steps/2, max.steps,1),]
      t13[t13!=0]=1; t13[t13==0]=0
      q_eta = mean(t13)
      
      coef_GXE_SS = c(coef_GXE_SS,q_eta)
    }
    
    
  }
  
  BI=max.steps/2
  coeff.E = apply(betas$GS.alpha[-(1:BI),,drop=FALSE], 2, stats::median);
  names(coeff.E) = E.names;
  coeff.clin = apply(betas$GS.b[-(1:BI),,drop=FALSE], 2, stats::median); 
  names(coeff.clin) = c(1,clin.names);
  
  coeff.G = apply(coef_G_RBLSS[-(1:BI),,drop=FALSE], 2, stats::median); 
  names(coeff.G) = G.names;
  coeff.GE = apply(coef_GXE_RBLSS[-(1:BI),], 2, stats::median); 
  names(coeff.GE) = GXE.names;
  
  coef_hat = c(coef_G_SS, coef_GXE_SS);names(coef_hat)=c(G.names,GXE.names)
  rank.GXE = sort(coef_hat,decreasing = TRUE)
  
  coefficient = list(clin=coeff.clin, E=coeff.E, GE=coeff.GE, G=coeff.G)
  out = list(GS.E=betas$GS.alpha, GS.C=betas$GS.b, GS.G=coef_G_RBLSS, GS.GXE=coef_GXE_RBLSS)
  fit = list(posterior = out, coefficient=coefficient, ranklist=rank.GXE, burn.in = BI, iterations=max.steps, design=list(g,xx=xx, CLC=c, E=e))
  
  fit
}
