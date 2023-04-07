
DataMatrix <- function(X, Y, E, clin,intercept=TRUE, debugging=FALSE)
{
  g = as.matrix(X);  y = Y
  n = nrow(g); p = ncol(g)
  c=clin; e=E
  #CLC=cbind(matrix(1,n,1),c)
  
  clin.names = E.names = G.names = NULL
  g = scale(g, center = TRUE, scale=FALSE)
  
  if(!is.null(y)){
    if(nrow(y) != n)  stop("Length of Y does not match the number of rows of X.");
  }
  
  if(!is.null(clin)){
    clin = as.matrix(clin)
    if(nrow(clin) != n)  stop("clin has a different number of rows than X.");
    if(is.null(colnames(clin))){colnames(clin)=paste("clin.", 1:ncol(clin), sep="")}
    CLC = clin
    noClin = FALSE
    clin.names = colnames(clin)
  }
  
  if(intercept){ # add intercept
    CLC = cbind(matrix(1,n,1,dimnames=list(NULL, "IC")), CLC)
  }
  
  if(!is.null(E)){
    E = as.matrix(E);env = ncol(E)
    E = scale(E, center = TRUE, scale=FALSE)
    if(nrow(E) != n)  stop("E has a different number of rows than X.");
    if(is.null(colnames(E))){colnames(E)=paste("E", 1:env, sep="")}
    E.names = colnames(E)
    #CLC = cbind(CLC, E)
    noE = FALSE
  }else if(!debugging){
    stop("E factors must be provided.")
  }
  
  if(is.null(colnames(g))){
    G.names = paste("G", 1:p, sep="")
  }else{
    G.names = colnames(g)
  }
  
  xx <- c()
  
  for (i in 1:p)
  {
    
    xx = cbind(xx,g[,i]*e)
    last = i*env; first = last-env+1
    #xx[,first:last] = cbind(x[,j], E*x[,j])
    colnames(xx)[first:last] = c(paste(G.names[i], "E", 1:env, sep=""))
    
  }
  
  xx <- scale(xx)
  
 
   
  
  dat = list(y=y, c=CLC, e=e, g=g, xx=xx, g=g,n=n, p=p, env=env, G.names=G.names, E.names=E.names, clin.names=clin.names, GXE.names=colnames(xx))
  return(dat) 
}
