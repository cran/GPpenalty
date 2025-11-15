
# log likelihood function
# default -----------------------
# no penalty
# profile log likelihood
# zero mean
# fixed nugget value with 1e-08
# no penalty
# separable covariance
# ------------------------------


nl <- function(par, y, x, lambda=0, d, mu=FALSE, g=FALSE, fixed_g=NULL, scad=FALSE, profile=TRUE){

  # kernel (need to define kernel here for foreach)
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


  eps <- sqrt(.Machine$double.eps)
  n <- length(y)
  theta <- par[1:d]
  if(profile) { # profile log-likelihood
    if(g) { # estimate g
      g_val <- par[d+1]
    } else { # fix g
      g_val <- if(!is.null(fixed_g)) fixed_g else eps
    }
    k <- kernel(x1=x, theta=theta, x2=NULL, g=g_val)
  } else { # regulr log-likelihood
    s2 <- par[d+1]
    if(mu) {
      if(g){
        mu_val <- par[d+2]
        g_val <- par[d+3]
      } else {
        mu_val <- par[d+2]
        g_val <- if(!is.null(fixed_g)) fixed_g else eps
      }
    } else {
      if(g) {
        mu_val <- 0
        g_val <- par[d+2]
      } else {
        mu_val <- 0
        g_val <- if(!is.null(fixed_g)) fixed_g else eps
      }
    }
    k <- s2*(kernel(x1=x, theta=theta, x2=NULL, g=g_val))
  }

  ki <- solve(k)

  if(profile) { # profile=T
    if (mu) { # solve for mu
      one <- rep(1, n)
      mu_val <- c(t(one) %*% ki %*% y / t(one) %*% ki %*% one)
    } else { # zero mean
      mu_val <- 0
    }
  }

  ldetk <- determinant(k, logarithm=TRUE)$modulus
  if(scad) { # scad penalty
    a <- 3.7
    scadtheta <- c()
    for(i in 1:length(theta)){
      if(theta[i] <= lambda) {
        scadtheta[i] <- lambda * theta[i]
      } else if(theta[i] > lambda & theta[i] <= a*lambda) {
        scadtheta[i] <- -(theta[i]^2-2*a*lambda*theta[i]+lambda^2)/(2*(a-1))
      } else if(theta[i] > a*lambda) {
        scadtheta[i] <- ((a+1)*lambda^2)/2
      }
    }
    if(profile) {
      ll <- -n/2*log((t(y-mu_val)%*%ki%*%(y-mu_val))) - (1/2)*(ldetk)-n*sum(scadtheta)
    } else {
      ll <- -(1/2)%*%(t(y-mu_val) %*% ki %*% (y-mu_val)) - (1/2)*(ldetk)-n*sum(scadtheta)
    }
  } else { # lasso penalty
    if(profile) {
      ll <- -n/2*log((t(y-mu_val)%*%ki%*%(y-mu_val))) - (1/2)*(ldetk)-n*lambda*sum(abs(theta))
    } else {
      ll <- -(1/2)%*%(t(y-mu_val) %*% ki %*% (y-mu_val)) - (1/2)*(ldetk)-n*lambda*sum(abs(theta))
    }
  }
  return(-ll)
}


# gradient
gradnl <- function(par, y, x, lambda=0, d, mu=FALSE, g=FALSE, fixed_g=NULL, scad=FALSE, profile=TRUE) {

  # kernel (need to define kernel here for foreach)
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


  # Euclidean distance (need to define this here for foreach)
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

  eps <- sqrt(.Machine$double.eps)
  n <- length(y)
  theta <- par[1:d]
  if(profile) { # profile
    if(g) { # estimate g
      g_val <- par[d+1]
    } else { # fix g
      g_val <- if(!is.null(fixed_g)) fixed_g else eps
    }
    k <- kernel(x1=x, theta=theta, x2=NULL, g=g_val)

  } else { # regular log-likelihood
    s2 <- par[d+1]
    if(mu) {
      if(g){
        mu_val <- par[d+2]
        g_val <- par[d+3]
      } else {
        mu_val <- par[d+2]
        g_val <- if(!is.null(fixed_g)) fixed_g else eps
      }
    } else {
      if(g) {
        mu_val <- 0
        g_val <- par[d+2]
      } else {
        mu_val <- 0
        g_val <- if(!is.null(fixed_g)) fixed_g else eps
      }
    }
    c <- kernel(x1=x, theta=theta, x2=NULL, g=g_val)
    k <- s2 * c
  }

  ki <- solve(k)
  one <- rep(1, n) # for mu

  if(profile) {
    if (mu) { # solve for mu
      mu_val <- c(t(one) %*% ki %*% y / t(one) %*% ki %*% one)
    } else { # zero mean
      mu_val <- 0
    }
  }

  dlltheta <- rep(NA, length(theta))

  kiy <- ki %*% (y-mu_val)

  #### theta ###
  # loop over theta
  for(i in 1:length(theta)) {
    if(profile) {
      dotk <- - k * euclidean_dist(x1=x[, i])
    } else {
      dotk <- -s2* c * euclidean_dist(x1=x[, i])
    }
    if(scad){ # scad
      a <- 3.7
      if(theta[i] <= lambda) {
        penalty <- lambda
      } else if(theta[i] > lambda & theta[i] <= a*lambda) {
        penalty <- -(theta[i]-a*lambda)/(a-1)
      } else if(theta[i] > a*lambda) {
        penalty <- 0
      }
      if(profile) {
        dlltheta[i] <- -1/2 * sum(diag(ki %*% dotk)) + n/2 * (t(kiy) %*% dotk %*% kiy) / (t(y-mu_val) %*% kiy) -
          n*penalty
      } else {
        dlltheta[i] <- -1/2 * sum(diag(ki %*% dotk)) + 1/2 * (t(kiy) %*% dotk %*% kiy) - n*penalty
      }


    } else { # lasso
      if(profile) {
        dlltheta[i] <- -1/2 * sum(diag(ki %*% dotk)) + n/2 * (t(kiy) %*% dotk %*% kiy) / (t(y-mu_val) %*% kiy) -
          n * lambda * sign(theta[i])
      } else {
        dlltheta[i] <- -1/2 * sum(diag(ki %*% dotk)) + 1/2 * (t(kiy) %*% dotk %*% kiy) -
          n * lambda * sign(theta[i])
      }
    }
  }

  # other parameters
  if(profile) {
    if(g) { # nugget
      dllg <- -1/2*sum(diag(ki))+ (n/2) * t(kiy) %*% kiy / (t(y-mu_val) %*% kiy)
      return(-c(dlltheta, dllg))
    } else {
      return(-c(dlltheta))
    }
  } else {
    # nugget
    dllg <- -1/2*s2*sum(diag(ki))+1/2* s2* t(kiy) %*% kiy
    # mu
    dllmu <- one %*% kiy
    # sigma^2
    dlls2 <- -1/2*sum(diag(ki %*% c)) + 1/2* t(kiy) %*% c %*% kiy

    # return
    if(mu) {
      if(g) {
        return(-c(dlltheta, dlls2, dllmu, dllg))
      } else {
        return(-c(dlltheta, dlls2, dllmu))
      }
    } else {
      if(g) {
        return(-c(dlltheta, dlls2, dllg))
      } else {
        return(-c(dlltheta, dlls2))
      }
    }
  }
}
