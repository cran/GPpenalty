# score function

#' @title score
#' @description Calculates a score value. Higher score values indicate better fits.
#'
#' @param y response variable vector
#' @param mu predicted mean vector
#' @param sigma predicted covariance matrix
#' @param md logical indicating whether to return to a Mahalanobis distance value (\code{md = TRUE}) and score value or
#'        only a score value (\code{md = FALSE})
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


score <- function(y, mu, sigma, md=FALSE) {
  R <- chol(sigma)
  logdet <- 2*sum(log(diag(R)))
  mh_dist <- crossprod(crossprod(solve(R), y - mu)) # mahalanobis distance
  mh_dist <- ifelse(mh_dist > 0 , mh_dist, 0) # if mahalanobis distance is negative, 0
  score <- (-logdet - mh_dist) / length(y)
  if(md) {
    return(list(score=c(score), md=c(mh_dist)))
  } else {
    return(c(score))
  }
}
