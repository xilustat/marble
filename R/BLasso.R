
BLasso <- function(X, Y, E, clin, max.steps, debugging=FALSE)
{
  dat = DataMatrix(X, Y, E, clin, intercept=TRUE, debugging=FALSE)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y;
  n = dat$n; p = dat$p; q=ncol(c)
  env = dat$env
  
  G.names = dat$G.names
  E.names = dat$E.names
  clin.names = dat$clin.names
  GXE.names = dat$GXE.names
  
  fun <- function(x)
  {
    pp = prod(x)
    if(sign(pp)==1) {1}
    else {0}
  }
  

hatAlpha = rep(1,env); hatb = rep(1,q); hatEta= rep(1,env); hatBeta=1
invSigAlpha0= diag(rep(1,env)); invSigb0 = diag(rep(1,q))
hatInvTauSq1= 1; hatInvTauSq2 = rep(1,env); hatLambdaSqStar1=1; hatLambdaSqStar2=1
hatSigmaSq=1; aStar=1; bStar=1; alpha=1; gamma=1
progress = 0

coef_G_BL = c(); coef_GXE_BL = c(); 

for(j in 1:p)
{
  
  x=g[,j]
  w=xx[,((env*(j-1)+1):(j*env))]
  betas =  BL(x,y,e,c,w,max.steps,n,hatBeta,hatEta,hatAlpha,hatb,hatInvTauSq1,hatInvTauSq2,invSigAlpha0, invSigb0, hatLambdaSqStar1, hatLambdaSqStar2,hatSigmaSq, aStar, bStar, alpha, gamma, progress)
  t1=as.matrix(betas$GS.beta)
  coef_G_BL=cbind(coef_G_BL, t1)
  
  for(d in 1:env)
  {
    t2 = as.matrix(betas$GS.eta[,d])
    coef_GXE_BL = cbind(coef_GXE_BL, t2)
  }
  
  
}


prob_seq=seq(0.2,1,by=0.05)
coef_G_SS_temp=c()

for (k in 1:length(prob_seq)) 
{
  beta_hat = c()
  prob=prob_seq[k]
  for (j in 1:ncol(coef_G_BL)) 
  {
    t11 = coef_G_BL[,j]
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
  for (j in 1:ncol(coef_GXE_BL)) 
  {
    t12 = coef_GXE_BL[,j]
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

coeff.G = apply(coef_G_BL[-(1:BI),,drop=FALSE], 2, stats::median); 
names(coeff.G) = G.names;
coeff.GE = apply(coef_GXE_BL[-(1:BI),], 2, stats::median); 
names(coeff.GE) = GXE.names;

coef_hat = c(coef_G_SS, coef_GXE_SS);names(coef_hat)=c(G.names,GXE.names)
rank.GXE = sort(coef_hat,decreasing = TRUE)

coefficient = list(clin=coeff.clin, E=coeff.E, GE=coeff.GE, G=coeff.G)
out = list(GS.E=betas$GS.alpha, GS.C=betas$GS.b, GS.G=coef_G_BL, GS.GXE=coef_GXE_BL)
fit = list(posterior = out, coefficient=coefficient, ranklist=rank.GXE, burn.in = BI, iterations=max.steps, design=list(g,xx=xx, CLC=c, E=e))

fit
}
