# score function

#' @title score
#' @description Calculates a score value. Higher score values indicate better fits.
#'
#' @param y response variable vector
#' @param mu predicted mean vector
#' @param sigma predicted covariance matrix
#' @param mh logical indicating whether to return to a Mahalanobis distance value (\code{mh = TRUE}) or
#'        a score value (\code{mh = FALSE})
#' @return a numerical value
#'
#'
#' @export
#'
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
#' ### get score value ###
#' score_value <- score(y.test, pred$mup, pred$Sigmap)
#'


score <- function(y, mu, sigma, mh=FALSE) {
  R <- chol(sigma)
  logdet <- 2*sum(log(diag(R)))
  mh_dist <- crossprod(crossprod(solve(R), y - mu)) # mahalanobis distance
  if(mh) {
    return(mh_dist)
  } else {
    mh_dist <- ifelse(mh_dist > 0 , mh_dist, 0) # if mahalanobis distance is negative, 0
    score <- (-logdet - mh_dist) / length(y)
    return(c(score))
  }
}
