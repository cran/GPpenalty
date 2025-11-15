# mle function

#' @title mle_gp
#' @description The function computes maximum likelihood estimates for the lengthscale, scale, mu, and nugget (g) parameters using \code{optim},
#'              with options to fix or assume zero for certain parameters.
#' @details The function uses numerical optimization for lengthscale and nugget parameters as
#'          there's no closed-form solutions. In contrast, closed form solutions exist for the scale and
#'          mu parameters. Users have options to choose whether to solve them analytically or include them in optimization process.
#'          If mu is assumed to be zero (by setting \code{mu=FALSE}), the input data should be centered beforehand.
#'          The nugget term (g) can also be optimized alongside the lengthscale parameter or fixed to a small constant.
#'          When no initial values are provided (\code{initialvals=NULL}), the function generates 10 random sets
#'          and selects the one that minimizes the negative log-likelihood. The number of sets can be specified by specifying \code{n_init}.
#'          Additionally, users can apply a penalty to the lengthscale parameter by specifying a tuning parameter, lambda.
#'          For guidance on choosing lambda, refer to \code{gp_cv} function.
#'
#' @param y A numeric vector of the response variable.
#' @param x A numeric vector or matrix of the input variables.
#' @param sep Logical indicator for using a separable kernel function (\code{sep=TRUE}) or an isotropic kernel function (\code{sep=FALSE}).
#'            Default is TRUE.
#' @param mu Logical indicator for assuming zero mean (\code{mu=FALSE}) or estimating the mean (\code{mu=TRUE}).
#'           Default is FALSE (assumes the data is centered beforehand).
#' @param g Logical indicator for fixing the nugget value to a small constant (\code{g=FALSE}) or estimating the nugget (\code{g=TRUE}). Default is FALSE.
#' @param fixed_g Nugget value to fix when \code{g=FALSE}. Default is \code{fixed_g=NULL}. If NULL, the nugget is fixed to 1.490116e-08.
#' @param profile Logical indicator for optimizing the profile log-likelihood (\code{profile=TRUE}). When TRUE, the log-likelihood is a function of lengthscale and nugget only.
#'                Solve the closed forms for scale and mu parameters. When FALSE, the full log-likelihood is optimized (lengthscale, scale, mean, and nugget are estimated together). Default is TRUE.
#' @param initialvals A numeric vector or matrix of initial values for optimization. The length should match the number of parameters to estimate.
#'                    Default is NULL. If NULL, 10 sets of initial values are randomly generated. The number of sets can be specified by specifying \code{n_init}.
#' @param n_init An integer indicating the number of randomly generated initial value sets to evaluate when \code{initialvals} is not provided.
#'               Default is 10.
#' @param penalty Logical indicator for penalization. Default is \code{penalty=FALSE} (returns MLE). When \code{penalty=TRUE} and no lambda value is specified, a set of estimated values along with evaluated lambda values is returned.
#' @param scad Logical indicator for a lasso penalty (\code{scad=FALSE}) or SCAD penalty (\code{scad=TRUE}) when \code{penalty=TRUE}. Default is lasso penalty.
#' @param lambda Tuning parameter value. Default is 0 (MLE). The user may specify a custom lambda value.
#' @param theta_upper Upper bound for theta in optim. Default is 1000.
#' @param theta_lower Lower bound for theta in optim. Default is 0.001.
#' @param ncores A number of cores for parallel computing with \code{optim}. Default is 1 (no parallelization). Make sure your system supports the specified number of cores.
#'
#' @return A list of y, x, and hyperparameters:
#' \itemize{
#'   \item \code{y}: A copy of y.
#'   \item \code{x}: A copy of x.
#'   \item \code{theta}: A matrix of estimated lengthscale parameter.
#'   \item \code{sigma2}: The estimated scale parameter.
#'   \item \code{mu}: Returns 0 if \code{mu=FALSE} otherwise the estimated mu parameter.
#'   \item \code{g}: Returns the \code{fixed_g} value if \code{g=FALSE} otherwise the estimated nugget value.
#'   \item \code{penalty}: A copy of the penalty indicator.
#'   \item \code{lambda}: A vector of evaluated lambda values if \code{penalty=TRUE} otherwise NULL.
#'   \item \code{theta_upper}: A copy of theta_upper for optimization.
#'   \item \code{theta_lower}: A copy of theta_lower for optimization.
#' }
#'
#'
#' @export
#'
#' @examples
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
#'
#' y <- f_x(x)
#'
#' ### Optimize only the lengthscale parameter and solve for scale. ###
#' ### Assume zero mean and fix g to a small constant. ###
#' fit <- mle_gp(y, x)
#'\donttest{
#' ### Include etimation of mu ###
#' fit <- mle_gp(y, x, mu=TRUE)
#'
#' ### Optimize g as well ###
#' fit <- mle_gp(y, x, mu=TRUE, g=TRUE)
#'
#' ### Jointly optimize the lengthscale and scale ###
#' fit <- mle_gp(y, x, profile=FALSE)
#'
#' ### Fix g to a user specified value ###
#' fit <- mle_gp(y, x, fixed_g=0.0001)
#'
#' ### Set the upper and lower bounds for theta ###
#' fit <- mle_gp(y, x, theta_upper=100, theta_lower=0.01)
#'
#' }
#'



# mle
mle_gp <- function(y, x, sep=TRUE, mu=FALSE, g=FALSE, fixed_g=NULL, profile=TRUE, initialvals=NULL, n_init = 10, penalty=FALSE, scad=FALSE, lambda=0, theta_upper=1000, theta_lower=0.001, ncores=1) {
  if(!is.null(dim(y))) y <- c(y) # y should be a vector
  if(is.list(x)) stop("x needs to be a vector or matrix")
  # check if x is a matrix
  if(is.null(dim(x))) x <- matrix(x, ncol=1)
  if(!is.null(dim(y))) y <- c(y) # y should be a vector
  # set dimension
  if(sep & ncol(x)!=1) {
    d <- ncol(x)
  } else {
    d <- 1
  }

  if(g & !is.null(fixed_g)) {
    fixed_g <- NULL
    warning("g will be estimated")
  }

  # if initial values are not specified
  if(is.null(initialvals)) {
    # # of parameters to estimate except for theta
    if(profile) {
      param <- if(g) 1 else 0 # if g=T then theta and g else only theta
    } else {
      if (mu & g) {
        param <- 3 # theta and (s2 + mu + g)
      } else if (mu | g) {
        param <- 2 # theta and (s2 + mu or g)
      } else {
        param <- 1 # theta and s2
      }
    }
    # generate initial values
    initialvals <- random_init(d, param, n_init)
  } else { # initial values are specified but the number doesn't match to the parameters to estimate
    if(is.null(dim(initialvals))) initialvals <- matrix(initialvals, ncol=1)
    if(profile & d!=ncol(initialvals)) {
      stop("the length of initial values doesn't match")
    }
    if(!profile) {
      n_init <- ncol(initialvals)
      if (mu & g) {
        param <- 3 # theta and (s2 + mu + g)
      } else if (mu | g) {
        param <- 2 # theta and (s2 + mu or g)
      } else {
        param <- 1 # theta and s2
      }
      if(n_init != d+param) {
        stop("the length of initial values doesn't match")
      }
    }
  }

  # if lambda value is specified, penalty is TRUE
  if(any(lambda!=0)) penalty <- TRUE

  object <- list()
  object$y <- y
  object$x <- x
  object$d <- d
  object$mu <- mu
  object$g <- g
  object$fixed_g <- fixed_g
  object$profile <- profile
  object$initialvals <- initialvals
  object$lambda <- lambda
  object$theta_upper <- theta_upper
  object$theta_lower <- theta_lower
  object$ncores <- ncores

  if(penalty) {
    estimates <- multi(object, penalty=penalty)
    return(list(y=y, x=x, theta=estimates$theta, sigma2=estimates$s2, mu=estimates$mu, g=estimates$g, penalty=penalty, lambda=estimates$lambda, theta_upper=theta_upper, theta_lower=theta_lower))
  } else {
    estimates <- multi(object)
    return(list(y=y, x=x, theta=estimates$theta, sigma2=estimates$s2, mu=estimates$mu, g=estimates$g, penalty=penalty, theta_upper=theta_upper, theta_lower=theta_lower))
  }

}
