
utils::globalVariables(c("g", "i"))

eps <- sqrt(.Machine$double.eps)



# kriging
kriging <- function(y, x, xx, theta, sigma2, mu=0, g=sqrt(.Machine$double.eps)) {

  Kx <- kernel(x1=x, theta=theta, x2=NULL, g=g)
  KXx <- kernel(x1=x, theta=theta, x2=xx, g=0)
  KXX <- kernel(x1=xx, theta=theta, x2=NULL, g=g)

  ki <- chol2inv(chol(Kx))
  mup <- mu + KXx %*% ki %*% (y-mu)
  R <- (KXX-crossprod(crossprod(solve(chol(Kx)),t(KXx))))
  Sigmap <- sigma2*R
  return(list(mup=mup, Sigmap=Sigmap, R=R))
}


# kernel
# kernel <- function(x1, theta, x2=NULL, nugget=NULL) {
#   # if(length(theta)==1 & ncol(x)!=1) theta <- rep(theta, ncol(x))
#   if(is.null(x2)) {
#     x2 <- x1
#     if(is.null(nugget)) nugget <- sqrt(.Machine$double.eps)
#     n <- nrow(x1)
#     # pairwise squared distance
#     k <- outer(1:n, 1:n, Vectorize(function(i, j) {
#       sum(theta*(x1[i,] - x2[j,])^2)
#     }))
#     k <- exp(-k) + diag(nugget, n)
#   } else {
#     n <- nrow(x1)
#     m <- nrow(x2)
#     k <- outer(1:n, 1:m, Vectorize(function(i, j) {
#       sum(theta*(x1[i,] - x2[j,])^2)
#     }))
#     k <- t(exp(-k))
#   }
#   return(k)
# }

# kernel function calls the c++ function
#' @title kernel
#' @description Compute the squared exponential kernel defined as \eqn{k = \exp(-\theta (x - x')^2) + g} , where \eqn{\theta} is the lengthscale parameter and \eqn{g} is a jitter term.
#'              Both isotropic and separable kernels are supported.
#' @details Matrix computations are implemented in C++ for improved performance and computational efficiency.
#'
#' @param x1 matrix of input locations
#' @param theta a scalar or vector specifying the lengthscale parameter. If a vector is provided, a separable kernel function is used.
#'              If a scalar is provided and \code{x1} has more than one column, an isotropic kernel is assumed.
#' @param x2 matrix of second input locations. If \code{NULL}, distance is computed between \code{x1} and itself.
#' @param g a jitter term. It is added when \code{x2=NULL} for computational stability.
#'
#' @return a matrix representing the evaluated kernel function
#' @useDynLib GPpenalty, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#'
#' @examples
#' ### isotropic ###
#' x <- matrix(seq(0, 10, length=10), ncol=1)
#' theta <- 5
#' k <- kernel(x1=x, theta=theta)
#'
#' ### anisotropic ###
#' x <- matrix(seq(0, 20, length=20), ncol=2)
#' theta <- c(2, 4)
#' k <- kernel(x1=x, theta=theta)
#'

kernel <- function(x1, theta, x2=NULL, g=NULL) {
  transpose <- F
  if(!is.matrix(x1)) x1 <- matrix(x1)

  if(!is.null(x2)) {
    transpose <- T
    if(!is.matrix(x2)) x2 <- matrix(x2)
    if(ncol(x1)!=ncol(x2)) {
      stop("the number of columns of x1 and x2 need to match")
    }
  }

  if(length(theta)==1 & ncol(x1)!=1) theta <- rep(theta, ncol(x1))
  if(is.null(x2)) {
    x2 <- x1
    if(is.null(g)) {
      g <- sqrt(.Machine$double.eps)
    }
  } else {
    if(is.null(g)) {
      g <- 0
    }
  }
  k <- kernel_exp(x1=x1, x2=x2, theta=theta, g=g)
  if(transpose) {
    k <- t(k)
  }
  return(k)
}

# Euclidean distance
# euclidean_dist <- function(x, X=NULL) {
#   x <- if(!is.matrix(x)) as.matrix(x)
#   if(is.null(X)) {
#     X <- x
#     n <- nrow(x)
#     # pairwise squared distance
#     k <- outer(1:n, 1:n, Vectorize(function(i, j) {
#       sum((x[i,] - X[j,])^2)
#     }))
#   } else {
#     X <- if(!is.matrix(X)) as.matrix(X)
#     n <- nrow(x)
#     m <- nrow(X)
#     k <- outer(1:n, 1:m, Vectorize(function(i, j) {
#       sum((x[i,] - X[j,])^2)
#     }))
#   }
#   return(k)
# }


# for c++ function
euclidean_dist <- function(x1, x2=NULL) {
  if(!is.matrix(x1)) x1 <- matrix(x1)

  if(!is.null(x2)) {
    if(!is.matrix(x2)) x2 <- matrix(x2)
    if(ncol(x1)!=ncol(x2)) {
      stop("the number of columns of x1 and x2 need to match")
    }
  }

  if(is.null(x2)) {
    x2 <- x1
  }
  k <- eucli_dist(x1=x1, x2=x2)

  return(k)
}

# generate initial values for optim randomly
#' @title random_init
#' @description generate random initial values for \code{optim}
#' @importFrom stats runif
#' @noRd
random_init <- function(d, n_param, n_init=10) {
  return(exp(matrix(runif(n_init*(d+n_param), -5, 2.5), nrow=n_init)))
}


# solve for scale parameter and mu variable
scale_mu <- function(y, x, theta, mu=FALSE, g=NULL) {
  # need to do like this for foreach
  kernel <- function(x1, theta, x2=NULL, g=NULL) {
    transpose <- F
    if(!is.matrix(x1)) x1 <- matrix(x1)

    if(!is.null(x2)) {
      transpose <- T
      if(!is.matrix(x2)) x2 <- matrix(x2)
      if(ncol(x1)!=ncol(x2)) {
        stop("the number of columns of x1 and x2 need to match")
      }
    }

    if(length(theta)==1 & ncol(x1)!=1) theta <- rep(theta, ncol(x1))
    if(is.null(x2)) {
      x2 <- x1
      if(is.null(g)) {
        g <- sqrt(.Machine$double.eps)
      }
    } else {
      if(is.null(g)) {
        g <- 0
      }
    }
    k <- kernel_exp(x1=x1, x2=x2, theta=theta, g=g)
    if(transpose) {
      k <- t(k)
    }
    return(k)
  }

  n <- length(y)
  if(is.null(g)) {
    g = eps
  }
  k <- kernel(x1=x, theta=theta, x2=NULL, g=g)
  ki <- solve(k)
  if(mu) {
    one <- rep(1, n)
    mu_val <- c(t(one) %*% ki %*% y / t(one) %*% ki %*% one)
    scale <- c((t(y-mu_val) %*% ki %*% (y-mu_val)) / n)
    return(list(scale=scale, mu=mu_val))
  } else {
    scale <- c((t(y) %*% ki %*% y) / n)
    return(list(scale = scale))
  }
}

# optimization function that evaluates 10 sets of initial values if users do not specify the initial values.

#' @title multi
#' @description Optimization function that evaluates 10 sets of initial values if users do not specify the initial values.
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster
#' @importFrom stats optim
#' @noRd

multi <- function(object, penalty=FALSE){


  if(is.null(object$lambda) & penalty) { # lambda is null and penalty=T
    lambda <- c(0, exp(seq(-7,2,length=40)))
  } else if(is.null(object$lambda) & !penalty) { # lambda is null and penalty=F
    lambda <- 0
  } else if(!(is.null(object$lambda)) & penalty) { # lambda is not null and penalty=T
    lambda <- object$lambda
  } else if(!(is.null(object$lambda)) & !penalty) { # lambda is not null and penalty=F
    lambda <- 0
  }
  if(!is.null(object$scad) & penalty) {
    scad <- object$scad
  } else { # default is lasso
    scad <- FALSE
  }
  if((!object$g & is.null(object$fixed_g)) | (is.null(object$fixed_g))) {
    fixed_g <- sqrt(.Machine$double.eps)
  } else if(object$fixed_g==-1) {
    fixed_g <- sqrt(.Machine$double.eps)
  } else {
    fixed_g <- object$fixed_g
  }

  y <- object$y
  x <- object$x
  d <- object$d
  mu <- object$mu
  g <- object$g
  profile <- object$profile
  theta_upper <- object$theta_upper
  theta_lower <- object$theta_lower

  initialvals <- as.matrix(object$initialvals)


  # upper bound and lower bound
  if(profile) {
    if(g) {
      upper <- c(rep(theta_upper, d), 1)
      lower <- c(rep(theta_lower, d), eps)
    } else {
      upper <- c(rep(theta_upper, d))
      lower <- c(rep(theta_lower, d))
    }
  } else {
    if(mu) {
      if(g) {
        upper <- c(rep(theta_upper, d), 100, 100, 1)
        lower <- c(rep(theta_lower, d), eps, -100, eps)
      } else {
        upper <- c(rep(theta_upper, d), 100, 100)
        lower <- c(rep(theta_lower, d), eps, -100)
      }
    } else {
      if(g) {
        upper <- c(rep(theta_upper, d), 100, 1)
        lower <- c(rep(theta_lower, d), eps, eps)
      } else {
        upper <- c(rep(theta_upper, d), 100)
        lower <- c(rep(theta_lower, d), eps)
      }
    }
  }
  # start parallel
  registerDoParallel(cores = object$ncores)

  estimates <- foreach(i=1:length(lambda),  .export = c("nl", "gradnl", "scale_mu")) %dopar% {
    lambda_val <- lambda[i]
    fit <-apply(initialvals,1,function(init){
      optim(init, fn=nl, gr=gradnl, method="L-BFGS-B",
            lower=lower, upper=upper, y=y, x=x, lambda=lambda_val, d=d, mu=mu, g=g, fixed_g=fixed_g, scad=scad, profile=profile)})
    # find the best
    min.ll<-which.min(unlist(lapply(fit,function(w){
      w[[2]] })))
    out <- fit[[min.ll]]

    # initialize
    theta <- rep(NA, d)
    ll <- NA
    conv <- NA
    s2 <- NA
    muval <- NA
    gval <- NA


    # parameters
    theta <- out$par[1:d]
    ll <- out$value
    conv <- out$convergence

    if(profile) {
      if(g) {
        gval <- out$par[d+1]
        param <- scale_mu(y, x, theta, mu=mu, g=gval)
        s2 <- param$scale
        if(mu) {
          muval <- param$mu
        } else {
          muval <- 0
        }
      } else {
        gval <- fixed_g
        param <- scale_mu(y, x, theta, mu=mu, g=gval)
        s2 <- param$scale
        if(mu) {
          muval <- param$mu
        } else {
          muval <- 0
        }
      }
    } else {
      s2 <- out$par[d+1]
      if(mu) {
        if(g) {
          muval <- out$par[d+2]
          gval <- out$par[d+3]
        } else {
          muval <- out$par[d+2]
          gval <- fixed_g
        }
      } else {
        muval <- 0
        if(g) {
          gval <- out$par[d+2]
        } else {
          gval <- fixed_g
        }
      }
    }
    list(theta=theta, ll=ll, conv=conv, s2=s2, muval=muval, gval=gval)

  }
  # stop parallel
  stopImplicitCluster()
  # parameters
  theta <- as.matrix(do.call(rbind, lapply(estimates, `[[`, "theta")), nrow=length(lambda))
  ll <- c(sapply(estimates, `[[`, "ll"))
  conv <- c(sapply(estimates, `[[`, "conv"))
  s2 <- c(sapply(estimates, `[[`, "s2"))
  muval <- c(sapply(estimates, `[[`, "muval"))
  gval <- c(sapply(estimates, `[[`, "gval"))

  if(penalty) {
    return(list(theta=theta, s2=s2, mu=muval, g=gval, ll=ll, conv=conv, lambda=lambda))
  } else {
    return(list(theta=theta, s2=s2, mu=muval, g=gval, ll=ll, conv=conv))
  }
}
