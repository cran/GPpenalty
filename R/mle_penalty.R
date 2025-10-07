# penalized mle

#' @title mle_penalty
#' @description Computes penalized maximum likelihood estimates for the lengthscale parameter using \code{optim}.
#' @details
#' This function takes the output from \code{\link{gp_cv}} and computes penalized MLEs for the lengthscale parameter,
#' along with MLEs for other model parameters. users may choose to apply the one standard error rule for selectingthe lambda value.
#' The \code{\link{gp_cv}} function returns both the optimal lambda and one standard error lambda. See \code{\link{gp_cv}} for details.
#'
#' @param object A list returned from \code{\link{gp_cv}}.
#' @param one.se Logical indicator for selecting the lambda value using the one-standard error. Default is FALSE. When FALSE, the lambda value that
#'               minimizes mse or maximizes score is selected. When TRUE, the lambda value is chosen based on the one-standard error rule.
#' @param lambda A user specified tuning parameter. This can be provided directly instead of performing cross-validation.
#' @param ncores A number of cores for parallel computing with \code{optim}. Default is 1 (no parallelization). Make sure your system supports the specified number of cores.
#'
#' @return A list of y, x, and hyperparameters:
#' \itemize{
#'   \item \code{y}: A copy of y.
#'   \item \code{x}: A copy of x.
#'   \item \code{theta}: A matrix of penalized lengthscale estimates.
#'   \item \code{sigma2}: The estimated scale parameter.
#'   \item \code{mu}: Returns 0 if \code{mu=FALSE} otherwise the estimated mu parameter.
#'   \item \code{g}: Returns the \code{fixed_g} value if \code{g=FALSE} otherwise the estimated nugget value.
#'   \item \code{lambda}: A scalar or vector of lambda values evaluated.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' ### training data ###
#' n <- 8
#'
#' ### test function ###
#' f_x <- function(x) {
#' return(sin(2*pi*x) + x^2)
#' }
#'
#' ### generate x ###
#' x <- runif(n, 0, 1)
#' y <- f_x(x)
#'
#' ### k-fold cross validation ###
#' cv.lambda <- gp_cv(y, x, k=4)
#'
#' ### fit the model ###
#' penalized.mle <- mle_penalty(cv.lambda)
#'
#' #### use the one standard error rule ###
#' penalized.mle <- mle_penalty(cv.lambda, one.se=TRUE)
#'
#' ### specify lambda ###
#' penalized.mle <- mle_penalty(cv.lambda, lambda=cv.lambda$lambda.score.max)
#'
#' }
#'

mle_penalty <- function(object, one.se=FALSE, lambda=NULL, ncores=1) {
  if(is.null(unname(unlist(object[grepl("lambda", names(object))]))) & is.null(lambda)) {
    stop("lambda value is not specified")
  } else if(is.null(lambda)) {
    if(object$metric=="score" & !one.se) {
      lambda <- unname(unlist(object[grepl("lambda.score.max", names(object))]))
    } else if(object$metric=="score" & one.se) {
      lambda <- unname(unlist(object[grepl("lambda.score.1se", names(object))]))
    } else if(object$metric=="mse" & !one.se) {
      lambda <- unname(unlist(object[grepl("lambda.min", names(object))]))
    } else if(object$metric=="mse" & one.se) {
      lambda <- unname(unlist(object[grepl("lambda.1se", names(object))]))
    }
  }

  if(is.null(dim(object$x))) object$x <- matrix(object$x, ncol=1)

  # if initial values are not specified
  if(is.null(object$initialvals)) {
    # # of parameters to estimate except for theta
    if(object$profile) {
      param <- if(g) 1 else 0
    } else {
      if (object$mu & object$g) {
        param <- 3 # theta and (s2 + mu + g)
      } else if (object$mu | object$g) {
        param <- 2 # theta and (s2 + mu or g)
      } else {
        param <- 1 # theta and s2
      }
    }
    # generate initial values
    initialvals <- random_init(object$dim, param)
  } else { # initial values are specified but the number doesn't match to the parameters to estimate
    if(is.null(object$initialvals)) object$initialvals <- matrix(object$initialvals, ncol=1)
    if(object$profile & object$dim!=ncol(object$initialvals)) {
      stop("the length of initial values doesn't match")
    }
    if(!object$profile) {
      n_init <- ncol(object$initialvals)
      if (object$mu & object$g) {
        param <- 3 # theta and (s2 + mu + g)
      } else if (object$mu | object$g) {
        param <- 2 # theta and (s2 + mu or g)
      } else {
        param <- 1 # theta and s2
      }
      if(n_init != object$dim+param) {
        stop("the length of initial values doesn't match")
      }
    }
  }

  object$lambda <- lambda
  object$ncores <- ncores

  estimates <- multi(object, penalty=TRUE)

  return(list(y=object$y, x=object$x, theta=estimates$theta, sigma2=estimates$s2, mu=estimates$mu, g=estimates$g, lambda=object$lambda))
}
