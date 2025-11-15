# decorrelated prediction error function

#' @title dpe
#' @description Calculates a decorrelated prediction error (DPE) value. Lower DPE values indicate better fits.
#'
#' @param y response variable vector
#' @param mu predicted mean vector
#' @param R predicted covariance matrix with the scale parameter removed
#' @return a numeric value
#'
#' @export
#' @examples
#'
#' ### test function ###
#' f_x <- function(x) {
#' return(sin(2*pi*x) + x^2)
#' }
#'
#' ### training data ###
#' n <- 8
#' x <- runif(n, 0, 1)
#' y <- f_x(x)
#'
#' ### testing data ###
#' n.test <- 100
#' x.test <- runif(n.test, 0, 1)
#' y.test <- f_x(x.test)
#'
#' ### get parameter estimates ###
#' out <- mle_gp(y, x)
#'
#' ### prediction ###
#' pred <- predict_gp(out, x.test)
#'
#' ### get DPE value ###
#' DPE_value <- dpe(y.test, pred$mup, pred$R)
#'



dpe <- function(y, mu, R) {
  chol_R <- chol(R)
  dpe <- crossprod(crossprod(solve(chol_R), y - mu))
  dpe <- ifelse(dpe>0, dpe, 0)
  return(c(sqrt(dpe)))
}
