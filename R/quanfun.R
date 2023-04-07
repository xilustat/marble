#' @importFrom stats quantile
quanfun <- function(x,prob)
{
  q_beta = as.matrix(quantile(x,c(prob,(1-prob))))
  pp = prod(q_beta)
  if(sign(pp)==1) {1}
  else {0}
}