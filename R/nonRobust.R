
nonRobust <- function(X, Y, E, clin, max.steps, sparse, debugging=FALSE)
{
  dat = DataMatrix(X, Y, E, clin, intercept=TRUE, debugging=FALSE)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y;
  n = dat$n; p = dat$p; q=ncol(c)
  env = dat$env
  
  #max.steps=10000
  hatAlpha = rep(1,env); hatb = rep(1,q); hatEta= rep(1,env); hatBeta=1
  invSigAlpha0= diag(rep(1,env)); invSigb0 = diag(rep(1,q))
  hatInvTauSq1= 1; hatInvTauSq2 = rep(1,env)
  sg1=1; sg2 = rep(1,env); hatLambdaSqStar1=1; hatLambdaSqStar2=1; hatSigmaSq=1
  aStar=1; bStar=1; alpha=1; gamma=1
  progress = 0; hatPiEta=1/2; hatPiBeta=1/2
  mu0=1; nu0=1
  
  
  if(sparse){
    result = BL_SS(X, Y, E, clin, max.steps, debugging=FALSE)
  }else{
    result = BLasso(X, Y, E, clin, max.steps, debugging=FALSE)
  }
  
  result
  
}