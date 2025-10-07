# cross validation function

#' @title gp_cv
#' @description Performs cross-validation to select an optimal tuning parameter for penalized MLE of the lengthscale parameter in Gaussian processes.
#' @details This function supports both leave-one-out and k-fold cross-validation for selecting a suitable tuning parameter value in penalized likelihoood estimation.
#'          Users can choose between evaluation metrics, including score and mse, to guide the selection process.
#'          The data is split into training and validation sets, and the model is trained on the training data and evaluated on both sets.
#'          This helps avoid selecting lambda values that lead to poor interpolation by the GP. The function returns the optimal lambda value along with the
#'          lambda selected using the one-standard error rule.
#'
#' @param y A numeric vector of the response variable.
#' @param x A numeric vector or matrix of the input variables.
#' @param lambda A tuning parameter. Default is NULL. Users may specify one or more lambda values to be evaluated.
#'               When NULL, 41 lambda values ranging from 0 to 7.389 will be automatically evaluated.
#' @param sep Logical indicator for using a separable kernel function (\code{sep=TRUE}) or an isotropic kernel function (\code{sep=FALSE}).
#'            Default is TRUE.
#' @param mu Logical indicator for assuming zero mean (\code{mu=FALSE}) or estimating the mean (\code{mu=TRUE}).
#'           Default is FALSE (assumes the data is centered beforehand).
#' @param g Logical indicator for fixing the nugget value to a small constant (\code{g=FALSE}) or estimating the nugget (\code{g=TRUE}). Default is FALSE.
#' @param fixed_g Nugget value to fix when \code{g=FALSE}. Default is \code{fixed_g=NULL}. If NULL, the nugget is fixed to 1.490116e-08.
#' @param profile Logical indicator for optimizing the profile log-likelihood (\code{profile=TRUE}). When TRUE, the log-likelihood is a function of lengthscale and nugget only.
#'                Solve the closed forms for scale and mu parameters. When FALSE, the full log-likelihood is optimized (lengthscale, scale, mean, and nugget are estimated together). Default is TRUE.
#' @param initialvals A numeric vector or matrix of initial values for optimization. The length should match the number of parameters to estimate.
#'                    Default is NULL. If NULL, 10 sets of initial values are randomly generated.
#' @param scad Logical indicator for a lasso penalty (\code{scad=FALSE}) or SCAD penalty (\code{scad=TRUE}) when \code{penalty=TRUE}. Default is lasso penalty.
#' @param k The number of folds for k-fold CV. Default is NULL. When NULL, leave-one-out CV using mean squared error metric is performed.
#'          To conduct k-fold CV, users must specify a value for \code{k}.
#' @param metric The evaluation metric used in CV. Default is \code{"score"}. The score metric is only available when \code{k} is specified.
#'               Supported metrics include score and mean squared error metrics. To use mean squared error metric, set \code{metric="mse"}.
#' @param ncores A number of cores for parallel computing with \code{optim}. Default is 1 (no parallelization). Make sure your system supports the specified number of cores.
#'               Paralleling is recommended to improve performance.
#'
#' @return A list includes y, x, selected lambda, and settings:
#' \itemize{
#'   \item \code{y}: A copy of y.
#'   \item \code{x}: A copy of x.
#'   \item \code{lambda.min}: Returned when \code{k} is not specified or \code{metric="mse"}; the lambda value that minimizes mean squared error across the folds.
#'   \item \code{lambda.1se}: Returned when \code{k} is not specified or \code{metric="mse"}; the lambda value selected using the one-standard-error rule.
#'   \item \code{lambda.score.max}: Returned when \code{k} is specified and \code{metric="score"}; the lambda value that maximizes the score across the folds.
#'   \item \code{lambda.score.1se}: Returned when \code{k} is specified and \code{metric="score"}; the lambda value selected using the one-standard-error rule.
#'   \item \code{initialvals}: A vector or matrix of initial values used in \code{optim}.
#'   \item \code{dim}: The dimensionality of the lengthscale parameter. If \code{sep=TRUE}, \code{dim} is equal to the number of columns in x. Otherwise it is set to 1 for isotropic kernels.
#'   \item \code{profile}: A copy of the logical indicator for profile likelihood optimization.
#'   \item \code{mu}: A copy of the logical indicator for mean estimation.
#'   \item \code{g}: A copy of the logical indicator for nugget estimation.
#'   \item \code{fixed_g}: The fixed nugget value used when \code{g = FALSE}. If NULL, the nugget is set to 1.490116e-08 in \code{mle_penalty} function.
#'   \item \code{metric}: A copy of the evaluation metric used in CV.
#'   \item \code{scad}: A copy of the logical indicator for SCAD penalty usage.
#' }
#' @useDynLib GPpenalty, .registration = TRUE
#' @importFrom Rcpp sourceCpp
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
#' }
#'


gp_cv <- function(y, x, lambda=NULL, sep=TRUE, mu=FALSE, g=FALSE, fixed_g=NULL, profile=TRUE, initialvals=NULL, scad=FALSE, k=NULL,
                             metric = "score", ncores=1) {
  if(!is.null(dim(y))) y <- c(y) # y should be a vector
  if(length(y)>30) message("This function is intended for small datasets and may take time to run.")
  if(is.list(x)) stop("x needs to be a vector or matrix")
  if(is.null(dim(x))) x <- matrix(x, ncol=1)
  if(sep & ncol(x)!=1) {
    dim <- ncol(x)
  } else {
    dim <- 1
  }
  n <- length(y)
  p <- ncol(x)
  # if k is not specified, leave one out
  if(is.null(k)) k <- length(y)
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
    initialvals <- random_init(dim, param)
  } else { # initial values are specified but the number doesn't match to the parameters to estimate
    if(profile & dim!=length(initialvals)) {
      stop("the length of initial values doesn't match")
    }
    if(!profile) {
      n_init <- length(initialvals)
      if (mu & g) {
        param <- 3 # theta and (s2 + mu + g)
      } else if (mu | g) {
        param <- 2 # theta and (s2 + mu or g)
      } else {
        param <- 1 # theta and s2
      }
      if(n_init != dim+param) {
        stop("the length of initial values doesn't match")
      }
    }
  }

  # tuning parameter
  if(is.null(lambda)) lambda <- c(0, exp(seq(-7,2,length=40)))


  # cv function in c++
  if(is.null(fixed_g)) {
    cv.lambda <- cv(y=y, x=x, lambda=lambda, n=n, k=k, p=p, dim=dim, initialvals=initialvals, profile=profile, mu=mu, g=g, scad=scad, ncores=ncores)
  } else {
    cv.lambda <- cv(y=y, x=x, lambda=lambda, n=n, k=k, p=p, dim=dim, initialvals=initialvals, profile=profile, mu=mu, g=g, fixed_g=fixed_g, scad=scad, ncores=ncores)
  }

  if(k==n) { # loocv
    metric = "mse"
    return(list(y=y, x=x, lambda.min = lambda[cv.lambda$mse_min], lambda.1se = lambda[cv.lambda$mse_1se], initialvals = initialvals, dim=dim, profile=profile, mu=mu, g=g, fixed_g=fixed_g, metric=metric, scad=scad))

  } else { # k-fold
    if(metric=="score") {
      return(list(y=y, x=x, lambda.score.max = lambda[cv.lambda$score_max], lambda.score.1se = lambda[cv.lambda$score_1se], initialvals = initialvals, dim=dim, profile=profile, mu=mu, g=g, fixed_g=fixed_g, metric=metric, scad=scad))
    } else {
      metric = "mse"
      return(list(y=y, x=x, lambda.min = lambda[cv.lambda$mse_min], lambda.1se = lambda[cv.lambda$mse_1se], initialvals = initialvals, dim=dim, profile=profile, mu=mu, g=g, fixed_g=fixed_g, metric=metric, scad=scad))
    }
  }
}

