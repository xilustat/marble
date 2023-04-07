
Robust <- function(X, Y, E, clin, max.steps, sparse, debugging)
{
  dat = DataMatrix(X, Y, E, clin, intercept=TRUE, debugging=FALSE)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y;
  n = dat$n; p = dat$p; q=ncol(c)
  env = dat$env
  
  hatAlpha = rep(1,env); hatb = rep(1,q); hatEta= rep(1,env); hatBeta=1; hatTau=1; hatV = rep(1,n)
  invSigAlpha0= diag(rep(1,env)); invSigb0 = diag(rep(1,q))
  hatSg1=1; hatSg2 = rep(1,env); hatEtaSq1=1; hatEtaSq2=1
  theta=0.5; r1=1; r=1; a=1; b=1; sg1=1; sg2 = rep(1,env)
  hatPiBeta=1/2; hatPiEta=1/2
  sh0=1; sh1=1
  progress = 0
  
  
  if(sparse){
    result = RobustBLSS(X, Y, E, clin, max.steps, debugging=FALSE)
  }else{
    result = RobustBL(X, Y, E, clin, max.steps, debugging=FALSE)
  }
  
  result
  
}
