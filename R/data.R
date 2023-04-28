#' simulated data for demonstrating the features of marble.
#'
#' Simulated gene expression data for demonstrating the features of marble.
#'
#' @docType data
#' @keywords datasets
#' @name dat
#' @aliases dat  X Y E clin
#' @usage data("dat")
#' @format dat consists of four components: X, Y, E, clin. 
#' @details
#'
#' \strong{The data model for generating Y}
#' 
#' Use subscript \eqn{i} to denote the \eqn{i}th subject. Let \eqn{(Y_{i}, X_{i}, E_{i}, clin_{i})} (\eqn{i=1,\ldots,n}) be
#' independent and identically distributed random vectors. \eqn{Y_{i}} is a continuous response variable representing the
#' phenotype. \eqn{X_{i}} is the \eqn{p}--dimensional vector of genetic factors. The environmental factors and clinical factors
#' are denoted as the \eqn{q}-dimensional vector \eqn{E_{i}} and the \eqn{m}-dimensional vector \eqn{clin_{i}}, respectively.
#' The \eqn{\epsilon} follows some heavy-tailed distribution. For \eqn{X_{ij}} (\eqn{j = 1,\ldots,p}), the measurement of the \eqn{j}th genetic factor on the \eqn{j}th subject, 
#' considering the following model:
#' \deqn{Y_{i} = \alpha_{0} + \sum_{k=1}^{q}\alpha_{k}E_{ik}+\sum_{t=1}^{m}\gamma_{t}clin_{it}+\beta_{j}X_{ij}+\sum_{k=1}^{q}\eta_{jk}X_{ij}E_{ik}+\epsilon_{i},}
#' where \eqn{\alpha_{0}} is the intercept, \eqn{\alpha_{k}}'s and \eqn{\gamma_{t}}'s are the regression coefficients corresponding to effects of environmental and clinical factors, respectively.
#' The \eqn{\beta_{j}}'s and \eqn{\eta_{jk}}'s are the regression coefficients of the genetic variants and G\eqn{\times}E interactions effects, correspondingly. 
#' The G\eqn{\times}E interactions effects are defined with \eqn{W_{j} = (X_{j}E_{1},\ldots,X_{j}E_{q}).} With a slight abuse of notation, denote \eqn{\tilde{W} = W_{j}.}  
#' Denote \eqn{\alpha=(\alpha_{1}, \ldots, \alpha_{q})^{T}}, \eqn{\gamma=(\gamma_{1}, \ldots, \gamma_{m})^{T}}, \eqn{\beta=(\beta_{1}, \ldots, \beta_{p})^{T}}, \eqn{\eta=(\eta_{1}^{T}, \ldots, \eta_{p}^{T})^{T}}, \eqn{\tilde{W} = (\tilde{W_{1}}, \dots, \tilde{W_{p}})}. 
#' Then model can be written as 
#' \deqn{Y_{i} = E_{i}\alpha + clin_{i}\gamma + X_{ij}\beta_{j} + \tilde{W}_{i}\eta_{j} + \epsilon_{i}.}
#'  
#' @examples
#' data(dat)
#' dim(X)
#' @seealso \code{\link{marble}}
NULL