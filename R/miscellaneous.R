
utils::globalVariables(c("g", "i"))

eps <- sqrt(.Machine$double.eps)



# kriging
kriging <- function(y, x, xx, theta, sigma2, mu=0, g=sqrt(.Machine$double.eps)) {

  Kx <- kernel(x1=x, theta=theta, nugget=g)
  KXx <- kernel(x1=x, x2=xx, theta=theta, nugget=0)
  KXX <- kernel(x1=xx, theta=theta, nugget=g)

  ki <- chol2inv(chol(Kx))
  mup <- mu + KXx %*% ki %*% (y-mu)
  Sigmap <- sigma2*(KXX-crossprod(crossprod(solve(chol(Kx)),t(KXx))))
  return(list(mup=mup, Sigmap=Sigmap))
}


# kernel
kernel <- function(x1, theta, x2=NULL, nugget=NULL) {
  # if(length(theta)==1 & ncol(x)!=1) theta <- rep(theta, ncol(x))
  if(is.null(x2)) {
    x2 <- x1
    if(is.null(nugget)) nugget <- sqrt(.Machine$double.eps)
    n <- nrow(x1)
    # pairwise squared distance
    k <- outer(1:n, 1:n, Vectorize(function(i, j) {
      sum(theta*(x1[i,] - x2[j,])^2)
    }))
    k <- exp(-k) + diag(nugget, n)
  } else {
    n <- nrow(x1)
    m <- nrow(x2)
    k <- outer(1:n, 1:m, Vectorize(function(i, j) {
      sum(theta*(x1[i,] - x2[j,])^2)
    }))
    k <- t(exp(-k))
  }
  return(k)
}

# Euclidean distance
euclidean_dist <- function(x, X=NULL) {
  x <- if(!is.matrix(x)) as.matrix(x)
  if(is.null(X)) {
    X <- x
    n <- nrow(x)
    # pairwise squared distance
    k <- outer(1:n, 1:n, Vectorize(function(i, j) {
      sum((x[i,] - X[j,])^2)
    }))
  } else {
    X <- if(!is.matrix(X)) as.matrix(X)
    n <- nrow(x)
    m <- nrow(X)
    k <- outer(1:n, 1:m, Vectorize(function(i, j) {
      sum((x[i,] - X[j,])^2)
    }))
  }
  return(k)
}


# generate initial values for optim randomly
#' @title random_init
#' @description generate random initial values for \code{optim}
#' @importFrom stats runif
#' @noRd
random_init <- function(dim, n_param) {
  return(exp(matrix(runif(10*(dim+n_param), -5, 2.5), nrow=10)))
}


# solve for scale parameter and mu variable
scale_mu <- function(y, x, theta, mu=FALSE, g=NULL) {
  # need to do like this for foreach
  kernel <- function(x1, theta, x2=NULL, nugget=NULL) {
    # if(length(theta)==1 & ncol(x)!=1) theta <- rep(theta, ncol(x))
    if(is.null(x2)) {
      x2 <- x1
      if(is.null(nugget)) nugget <- sqrt(.Machine$double.eps)
      n <- nrow(x1)
      # pairwise squared distance
      k <- outer(1:n, 1:n, Vectorize(function(i, j) {
        sum(theta*(x1[i,] - x2[j,])^2)
      }))
      k <- exp(-k) + diag(nugget, n)
    } else {
      n <- nrow(x1)
      m <- nrow(x2)
      k <- outer(1:n, 1:m, Vectorize(function(i, j) {
        sum(theta*(x1[i,] - x2[j,])^2)
      }))
      k <- t(exp(-k))
    }
    return(k)
  }

  n <- length(y)
  if(is.null(g)) {
    g = eps
  }
  k <- kernel(x1=x, theta=theta, nugget=g)
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
  dim <- object$dim
  mu <- object$mu
  g <- object$g
  profile <- object$profile
  initialvals <- object$initialvals



  # upper bound and lower bound
  if(profile) {
    if(g) {
      upper <- c(rep(1000, dim), 1)
      lower <- c(rep(1/1000, dim), eps)
    } else {
      upper <- c(rep(1000, dim))
      lower <- c(rep(1/1000, dim))
    }
  } else {
    if(mu) {
      if(g) {
        upper <- c(rep(1000, dim), 100, 100, 1)
        lower <- c(rep(1/1000, dim), eps, -100, eps)
      } else {
        upper <- c(rep(1000, dim), 100, 100)
        lower <- c(rep(1/1000, dim), eps, -100)
      }
    } else {
      if(g) {
        upper <- c(rep(1000, dim), 100, 1)
        lower <- c(rep(1/1000, dim), eps, eps)
      } else {
        upper <- c(rep(1000, dim), 100)
        lower <- c(rep(1/1000, dim), eps)
      }
    }
  }
  # start parallel
  registerDoParallel(cores = object$ncores)

  estimates <- foreach(i=1:length(lambda), .export = c("nl", "gradnl", "scale_mu")) %dopar% {
    lambda_val <- lambda[i]
    fit <-apply(initialvals,1,function(init){
      optim(init, fn=nl, gr=gradnl, method="L-BFGS-B",
            lower=lower, upper=upper, y=y, x=x, lambda=lambda_val, dim=dim, mu=mu, g=g, fixed_g=fixed_g, scad=scad, profile=profile)})
    # find the best
    min.ll<-which.min(unlist(lapply(fit,function(w){
      w[[2]] })))
    out <- fit[[min.ll]]

    # initialize
    theta <- rep(NA, dim)
    ll <- NA
    conv <- NA
    s2 <- NA
    muval <- NA
    gval <- NA


    # parameters
    theta <- out$par[1:dim]
    ll <- out$value
    conv <- out$convergence

    if(profile) {
      if(g) {
        gval <- out$par[dim+1]
        param <- scale_mu(y, x, theta, mu=mu, g=gval)
        s2 <- param$scale
        if(mu) {
          muval <- param$mu
        } else {
          muval <- 0
        }
        # initialvals[1, ] <- c(theta[i, ], gval[i])
      } else {
        gval <- fixed_g
        param <- scale_mu(y, x, theta, mu=mu, g=gval)
        s2 <- param$scale
        if(mu) {
          muval <- param$mu
        } else {
          muval <- 0
        }
        # initialvals[1, ] <- c(theta[i, ])
      }
    } else {
      s2 <- out$par[dim+1]
      if(mu) {
        if(g) {
          muval <- out$par[dim+2]
          gval <- out$par[dim+3]
          # initialvals[1, ] <- c(theta[i, ], s2[i], muval[i], gval[i])

        } else {
          muval <- out$par[dim+2]
          # initialvals[1, ] <- c(theta[i, ], s2[i], muval[i])
          gval <- fixed_g
        }
      } else {
        muval <- 0
        if(g) {
          gval <- out$par[dim+2]
          # initialvals[1, ] <- c(theta[i, ], s2[i], gval[i])
        } else {
          # initialvals[1, ] <- c(theta[i, ], s2[i])
          gval <- fixed_g
        }
      }
    }
    list(theta=theta, ll=ll, conv=conv, s2=s2, muval=muval, gval=gval)

  }
  # stop parallel
  stopImplicitCluster()
  # parameters
  theta <- do.call(rbind, lapply(estimates, `[[`, "theta"))
  ll <- sapply(estimates, `[[`, "ll")
  conv <- sapply(estimates, `[[`, "conv")
  s2 <- sapply(estimates, `[[`, "s2")
  muval <- sapply(estimates, `[[`, "muval")
  gval <- sapply(estimates, `[[`, "gval")

  if(penalty) {
    return(list(theta=theta, s2=s2, mu=muval, g=gval, ll=ll, conv=conv, lambda=lambda))
  } else {
    return(list(theta=theta, s2=s2, mu=muval, g=gval, ll=ll, conv=conv))
  }
}
