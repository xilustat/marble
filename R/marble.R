#' @useDynLib marble, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#' fit a robust Bayesian variable selection model for G×E interactions.
#' @keywords models
#' @param X the matrix of predictors (genetic factors). Each row should be an observation vector. 
#' @param Y the continuous response variable. 
#' @param E a matrix of environmental factors. E will be centered. The interaction terms between X (genetic factors) and E will be automatically created and included in the model.
#' @param clin a matrix of clinical variables. Clinical variables are not subject to penalty. Clinical variables will be centered and a column of 1 will be added to the Clinical matrix as the intercept.
#' @param max.steps the number of MCMC iterations.
#' @param robust logical flag. If TRUE, robust methods will be used.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly.
#' @param debugging logical flag. If TRUE, progress will be output to the console and extra information will be returned.
#' @return an object of class `marble' is returned, which is a list with component:
#' \item{posterior}{the posterior samples of coefficients from the MCMC.}
#' \item{coefficient}{the estimated value of coefficients.}
#' \item{ranklist}{the rank list of main and interaction effects.}
#' \item{burn.in}{the total number of burn-ins.}
#' \item{iterations}{the total number of iterations.}
#' \item{design}{the design matrix of all effects.}
#' 
#' @details Consider the data model described in "\code{\link{dat}}":
#' \deqn{Y_{i} = \alpha_{0} + \sum_{k=1}^{q}\alpha_{k}E_{ik}+\sum_{t=1}^{m}\gamma_{t}clin_{it}+\beta_{j}X_{ij}+\sum_{k=1}^{q}\eta_{jk}X_{ij}E_{ik}+\epsilon_{i},}
#' Where \eqn{\alpha_{0}} is the intercept, \eqn{\alpha_{k}}'s and \eqn{\gamma_{t}}'s are the regression coefficients corresponding to effects of environmental and clinical factors.
#' And \eqn{\beta_{j}}'s and \eqn{\eta_{jk}}'s are the regression coefficients of the genetic variants and G\eqn{\times}E interactions effects, correspondingly. 
#'
#' When sparse=TRUE (default), spike--and--slab priors are imposed to identify important main and interaction effects. If sparse=FALSE, Laplacian shrinkage will be used.
#'
#' When robust=TRUE (default), the distribution of \eqn{\epsilon_{i}} is defined as a Laplace distribution with density
#' \eqn{
#' f(\epsilon_{i}|\nu) = \frac{\nu}{2}\exp\left\{-\nu |\epsilon_{i}|\right\}
#' }, (\eqn{i=1,\dots,n}), which leads to a Bayesian formulation of LAD regression. If robust=FALSE, \eqn{\epsilon_{i}} follows a normal distribution.
#'
#' Here, a rank list of the main and interaction effects is provided. For method incorporating spike-and-slab priors,  
#' the inclusion probability is used to indicate the importance of predictors. 
#' We use a binary indicator \eqn{\phi} to denote that the membership of the non-spike distribution. 
#' Take the main effect of the \eqn{j}th genetic factor, \eqn{X_{j}}, as an example. 
#' Suppose we have collected H posterior samples from MCMC after burn-ins. The \eqn{j}th G factor is included 
#' in the marginal G\eqn{\times}E model at the \eqn{j}th MCMC iteration if the corresponding indicator is 1, i.e., \eqn{\phi_j^{(h)} = 1}. 
#' Subsequently, the posterior probability of retaining the \eqn{j}th genetic main effect in the final marginal model is defined as the average of all the indicators for the \eqn{j}th G factor among the H posterior samples. 
#' That is, \eqn{p_j = \hat{\pi} (\phi_j = 1|y) = \frac{1}{H} \sum_{h=1}^{H} \phi_j^{(h)}, \;  j = 1, \dots,p.}
#' A larger posterior inclusion probability \eqn{j}th indicates a stronger empirical evidence that the \eqn{j}th genetic main effect has a non-zero coefficient, i.e., a stronger association with the phenotypic trait. 
#' For method without spike-and-slab priors, variable selection is based on different level of credible intervals.
#' 
#' Both \eqn{X}, \eqn{clin} and \eqn{E} will be standardized before the generation of interaction terms to avoid the multicollinearity between main effects and interaction terms.
#'
#' Please check the references for more details about the prior distributions.
#' 
#' @references
#' Lu, X., Fan, K., Ren, J., and Wu, C. (2021). Identifying Gene–Environment Interactions With Robust Marginal Bayesian Variable Selection.
#' {\emph{Frontiers in Genetics}, 12:667074} \doi{10.3389/fgene.2021.667074}
#'
#' @seealso \code{\link{GxESelection}}
#'
#' @examples
#' data(dat)
#'
#' ## default method
#' max.steps=5000
#' fit=marble(X, Y, E, clin, max.steps=max.steps)
#' 
#' ## coefficients of parameters
#' fit$coefficient
#'
#' ## Estimated values of main G effects 
#' fit$coefficient$G
#' 
#' ## Estimated values of interactions effects 
#' fit$coefficient$GE
#'
#' ## Rank list of main G effects and interactions 
#' fit$ranklist
#'
#' \donttest{
#' ## alternative: robust selection
#' fit=marble(X, Y, E, clin, max.steps=max.steps, robust=TRUE, sparse=FALSE)
#' fit$coefficient
#' fit$ranklist
#'
#' ## alternative: non-robust sparse selection
#' fit=marble(X, Y, E, clin, max.steps=max.steps, robust=FALSE, sparse=FALSE)
#' fit$coefficient
#' fit$ranklist
#' }
#'
#' @export
marble <- function(X, Y, E, clin, max.steps=10000, robust=TRUE, sparse=TRUE, debugging=FALSE)
{

  dat = DataMatrix(X, Y, E, clin, intercept=TRUE, debugging=FALSE)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y;
  n = dat$n; p = dat$p; q=ncol(c)
  env = dat$env
  
  G.names = dat$G.names
  E.names = dat$E.names
  clin.names = dat$clin.names
  GXE.names = dat$GXE.names
  
  
  if(robust){
    out = Robust(X, Y, E, clin, max.steps, sparse, debugging)
  }else{
    out = nonRobust(X, Y, E, clin, max.steps, sparse, debugging)
  }
  class(out) = "marble"
  out
  
}
