#' print a marble object
#' 
#' Print a summary of a marble object
#' 
#' @param x marble object.
#' @param digits significant digits in printout.
#' @param ... other print arguments.
#' @return No return value, called for side effects.
#' @usage \method{print}{marble}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{marble}} 
#' @export
print.marble=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCoefficients:\n")
  print(x$coefficient, digits)
}
#'
#' print a GxESelection object
#' 
#' Print a summary of a GxESelection object
#' 
#' @param x GxESelection object.
#' @param digits significant digits in printout.
#' @param ... other print arguments.
#' @return No return value, called for side effects.
#' @usage \method{print}{GxESelection}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{GxESelection}}
#' @export
print.GxESelection=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nMethod:\n")
  print(x$method)
  cat("\n")
  print(x$effects)
}
