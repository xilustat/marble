#' Variable selection for a marble object
#' 
#' @param obj marble object.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly.
#' @details For class `Sparse',  the inclusion probability is used to indicate the importance of predictors. 
#' Here we use a binary indicator \eqn{\phi} to denote that the membership of the non-spike distribution. 
#' Take the main effect of the \eqn{j}th genetic factor, \eqn{X_{j}}, as an example. 
#' Suppose we have collected H posterior samples from MCMC after burn-ins. The \eqn{j}th G factor is included 
#' in the marginal G\eqn{\times}E model at the \eqn{j}th MCMC iteration if the corresponding indicator is 1, i.e., \eqn{\phi_j^{(h)} = 1}. 
#' Subsequently, the posterior probability of retaining the \eqn{j}th genetic main effect in the final marginal model is defined as the average of all the indicators for the \eqn{j}th G factor among the H posterior samples. 
#' That is, \eqn{p_j = \hat{\pi} (\phi_j = 1|y) = \frac{1}{H} \sum_{h=1}^{H} \phi_j^{(h)}, \;  j = 1, \dots,p.}
#' A larger posterior inclusion probability of \eqn{j}th indicates a stronger empirical evidence that the \eqn{j}th genetic main effect has a non-zero coefficient, i.e., a stronger association with the phenotypic trait. 
#' Here, we use 0.5 as a cutting-off point. If \eqn{p_j > 0.5}, then the \eqn{j}th genetic main effect is included in the final model. Otherwise, the \eqn{j}th genetic main effect is excluded in the final model.
#' For class `NonSparse', variable selection is based on 95\% credible interval.
#' Please check the references for more details about the variable selection.
#' 
#' @references
#' Lu, X., Fan, K., Ren, J., and Wu, C. (2021). Identifying Gene–Environment Interactions With Robust Marginal Bayesian Variable Selection.
#' {\emph{Frontiers in Genetics}, 12:667074} \doi{10.3389/fgene.2021.667074}
#' 
#' @rdname GxESelection
#' @return an object of class `GxESelection' is returned, which is a list with components:
#' \item{method}{method used for identifying important effects.}
#' \item{effects}{a list of indicators of selected effects.}
#' 
#' @seealso \code{\link{marble}}
#' @examples
#' data(dat)
#' max.steps=5000
#' ## sparse
#' fit=marble(X, Y, E, clin, max.steps=max.steps)
#' selected=GxESelection(fit,sparse=TRUE)
#' selected
#'
#' ## non-sparse
#' fit=marble(X, Y, E, clin, max.steps=max.steps, sparse=FALSE)
#' selected=GxESelection(fit,sparse=FALSE)
#' selected
#' 
#'
#' 
#' @export
GxESelection <- function(obj,sparse){
  if(sparse){
    out = GxESelection.Sparse(obj)
  }else{
    out = GxESelection.NonSparse(obj, prob=0.95)
  }
  
  out
}

