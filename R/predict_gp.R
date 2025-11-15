# prediction

#' @title predict_gp
#' @description Computes the posterior mean and covariance matrix for a given set of
#'              input locations based on a fitted model.
#' @details
#' From the model fitted by \code{\link{mle_gp}} or \code{\link{mle_gp}}, the posterior mean and
#' covariance matrix are computed.
#'
#' @param out out from \code{\link{mle_gp}} or \code{\link{mle_gp}}.
#' @param xx A numerical vector or matrix of new input locations.
#'
#' @return A list of predictive posterior mean and covariance:
#' \itemize{
#'   \item \code{mup}: vector of predicted posterior mean
#'   \item \code{Sigmap}: predictive posterior covariance matrix
#'   \item \code{R}: predictive posterior covariance matrix with the scale parameter removed
#' }
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
#'
predict_gp <- function(out, xx) {
  if(is.list(xx)) stop("xx needs to be a vector or matrix")
  if(is.null(dim(xx))) xx <- matrix(xx, ncol=1)
  # d <- out$dim
  y <- out$y
  x <- if(!is.null(dim(out$x))) out$x else matrix(out$x, ncol=1)
  theta <- out$theta
  mu <- if(!is.null(out$mu)) out$mu else 0
  g <- if(!is.null(out$g)) out$g else sqrt(.Machine$double.eps)
  if(is.null(out$sigma2)) {
    sigma2 <- scale_mu(y, x, theta)
  } else {
    sigma2 <- out$sigma2
  }
  results <- kriging(y, x, xx, theta, sigma2, mu, g)
  return(results)
}



