# this file contains functions that probably do not need to be directly accessed
# but are necessary for the functions in MZIP.R to work. 


theta.to.params <- function(theta, n, adj_mat=0, marginal=F, arp=1, parp=1){
  
  # helper function for mzip.sim and mzip.fitted
  # returns a list of vectors and matrices of parameters separated back into d, A, B, gamma, alpha, beta
  
  adj <- any(adj_mat!=0)
  nmat <- sum(adj_mat)
  d <- theta[1:n]
  if(adj==T){
    A <- restrict_mat(theta[(n + 1):(n + nmat)], adj_mat)
    B <- restrict_mat(theta[(n + nmat + 1):(n + (arp+1)*nmat)], adj_mat)
    if (marginal){
      pd <- theta[(n + (arp+1)*nmat + 1): (2*n + (arp+1)*nmat)]
      pa <- restrict_mat(theta[(2*n + (arp+1)*nmat + 1):(2*n + (arp+2)*nmat)], adj_mat)
      pb <- restrict_mat(theta[(2*n + nmat*(arp+2)+ 1): (2*n + nmat*(arp+parp+2))], adj_mat)
    } else{
      pd <- theta[n + (arp+1)*nmat + 1]
      pa <- theta[(n + (arp+1)*nmat + 2): (2*n + (arp+1)*nmat + 1)]
      pb <- theta[(2*n + (arp+1)*nmat + 2): ((2+parp)*n + (arp+1)*nmat + 1)]
    }
  } else{
    A <- matrix(theta[(n + 1):(n + n^2)], nrow=n, byrow=TRUE)
    B <- matrix(theta[(n + n^2 + 1):(n + (arp+1)*n^2)], nrow=n, byrow=TRUE)
    if (marginal){
      pd <- theta[(n + (arp+1)*n^2 + 1) : (2*n + (arp+1)*n^2)]
      pa <- matrix(theta[(2*n + (arp+1)*n^2 + 1):(2*n + (arp+parp+1)*n^2)], nrow=n, byrow=T)
      pb <- matrix(theta[(2*n + (arp+parp+1)*n^2 + 1): (2*n + (arp+parp+2)*n^2)], nrow=n, byrow=T)
    } else{
      pd <- theta[n + (arp+1)*n^2 + 1]
      pa <- theta[(n + (arp+1)*n^2 + 2):(2*n + (arp+1)*n^2 + 1)]
      pb <- theta[(2*n + (arp+1)*n^2 + 2): ((2+parp)*n + (arp+1)*n^2 + 1)]
    }
  }
  return (list("d"=d, "A"=A, "B"=B, "pd"=pd, "pa"=pa, "pb"=pb))
}


restrict_mat <- function(values, adj_mat) {
  
  # helper function for imposing adjacency matrix constraint
  
  if (sum(adj_mat) != length(values)) {
    stop("Number of values does not match the count of 1s in the adjacency matrix.")}
  n <- nrow(adj_mat)
  mat <- matrix(0, nrow = n, ncol = n)
  coordinates <- which(adj_mat == 1, arr.ind = TRUE)
  for (i in 1:nrow(coordinates)) {
    mat[coordinates[i, 1], coordinates[i, 2]] <- values[i]
  }
  return(mat)
}

fill_diag <- function(vals, adj_mat){
  
  # helper function for initialising values in mzip.EM when the parameter space is
  # restricted to where there are 1's in the adjacency matrix
  
  n <- length(adj_mat[1, ])
  vec <- c()
  for (i in 1:n){
    for (j in 1:n){
      if (adj_mat[i, j] == 1){
        if (i==j){
          vec <- c(vec, vals[i])
        } else {
          vec <- c(vec, 0)
        }
      }
    }
  }
  return(vec)
}


mzip.fitted.o <- function(theta, data, ll=F, marginal=F, arp=1, parp=1, adj_mat=0){
  
  # Part of mzip.fitted.
  # Returns a list of fitted intensities and probabilities of zero-inflation,
  # based on the estimated parameters [theta] and the previous data.
  # The intensities and probabilities at time 1 are initialised at the intercept terms d and gamma.
  
  n <- length(data[, 1])
  t.len <- length(data[1, ])
  maxar <- max(arp, parp)
  
  z <- theta.to.params(theta, n, marginal=marginal, arp=arp, parp=parp, adj_mat=adj_mat)
  d <- z$d; A <- z$A; B <- z$B; pd <- z$pd; pa <- z$pa; pb <- z$pb
  
  # set up vector of relevant past values for t = maxar + 1
  count.prev <- c()
  for (t in 1:arp){count.prev <- c(data[, t], count.prev)}
  p.count.prev <- c()
  for (t in 1:parp){p.count.prev <- c(data[, t], p.count.prev)}
  
  intensities <- matrix(0, nrow=n, ncol=t.len)
  if (marginal){
    probs <- matrix(0, nrow=n, ncol=t.len)
    probs[, 1:maxar] <- 1 / (1 + exp(-pd))
  } else{
    probs <- rep(0, t.len)
    probs[1:maxar] <- 1 / (1 + exp(-pd))
  }
  intensities[, 1:maxar] <- d
  
  for (t in (maxar+1):t.len){
    if (ll==T){
      intensities[, t] <- d + A %*% intensities[, t-1] + B %*% log(count.prev + 1)
      if (marginal){
        probs[, t] <- 1 / (1 + exp(-(pd + pa %*% exp(intensities[, t-1]) + pb %*% p.count.prev)))
      } else{
        probs[t] <- 1 / (1 + exp(-(pd + sum(pa*exp(intensities[, t-1])) + sum(pb*p.count.prev))))
      }
    }
    else {
      intensities[, t] <- d + A %*% intensities[, t-1] + B %*% count.prev
      if (marginal){
        probs[, t] <- 1 / (1 + exp(-(pd + pa %*% intensities[, t-1] + pb %*% p.count.prev)))
      } else{
        probs[t] <- 1 / (1 + exp(-(pd + sum(pa*intensities[, t-1]) + sum(pb*p.count.prev))))
      }
    }
    
    if (arp==1) {
      count.prev <- data[, t]
    } else{
      count.prev <- c(data[, t], count.prev[1:((arp-1)*n)])
    }
    if (parp==1) {
      p.count.prev <- data[, t]
    } else{
      p.count.prev <- c(data[, t], p.count.prev[1:((parp-1)*n)])
    }
  }
  if (ll==T){
    intensities <- exp(intensities)
  }
  return (list(ints=intensities, probs=probs))
}


mzip.fitted.g <- function(theta, data, adj_mat, ll=F, exo=0, global.neighbour=T,
                          global.self=F){
  # Part of mzip.fitted.
  # Returns a list of fitted intensities and probabilities of zero-inflation,
  # for the more constrained marginally zero-inflated model. 
  
  n <- length(data[, 1])
  t.len <- length(data[1, ])
  if (global.neighbour){
    d <- theta[1:n]; a1 <- theta[(n+1):(2*n)]; a2 <- theta[2*n+1];
    b1 <- theta[(2*n+2):(3*n+1)]; b2 <- theta[3*n+2]
    pd <- theta[(3*n+3):(4*n+2)]; pa1 <- theta[(4*n+3):(5*n+2)]; pa2 <- theta[5*n+3]
    pb1 <- theta[(5*n+4):(6*n+3)]; pb2 <- theta[6*n+4]
    if (any(exo != 0)){
      c <- theta[(6*n+5):(7*n+4)]
    }
    
  } else{
    z <- matrix(theta, nrow=n)
    d <- z[, 1]; a1 <- z[, 2]; a2 <- z[, 3]; b1 <- z[, 4]; b2 <- z[, 5]
    pd <- z[, 6]; pa1 <- z[, 7]; pa2 <- z[, 8]; pb1 <- z[, 9]; pb2 <- z[, 10]
    if (any(exo != 0)){
      c <- z[, 11]
    }
  }
  
  ints <- matrix(0, n, t.len)
  ints[, 1] <- d
  probs <- matrix(0, n, t.len)
  probs[, 1] <- 1 / (1 + exp(-pd))
  
  if (ll){
    for (t in 2:t.len){
      if (any(exo != 0)){
        ints[, t] <- d + a1*ints[, t-1] + a2*adj_mat%*%ints[, t-1] + b1*log(data[, t-1]+1) + b2*adj_mat%*%log(data[, t-1]+1) + c*exo[, t-1]
      } else{
        ints[, t] <- d + a1*ints[, t-1] + a2*adj_mat%*%ints[, t-1] + b1*log(data[, t-1]+1) + b2*adj_mat%*%log(data[, t-1]+1)
      }
      probs[, t] <- 1 / (1 + exp(-(pd + pa1*exp(ints[, t-1]) + pa2*adj_mat%*%exp(ints[, t-1]) + pb1*data[, t-1] + pb2*adj_mat%*%data[, t-1])))
    }
    ints <- exp(ints)
  }
  else {
    for (t in 2:t.len){
      ints[, t] <- d + a1*ints[, t-1] + a2*adj_mat%*%ints[, t-1] + b1*data[, t-1] + b2*adj_mat%*%data[, t-1]
      probs[, t] <- 1 / (1 + exp(-(pd + pa1*ints[, t-1] + pa2*adj_mat%*%ints[, t-1] + pb1*data[, t-1] + pb2*adj_mat%*%data[, t-1])))
    }
  }
  return(list('ints'=ints, 'probs'=probs))
}


mzip.Mstep <- function(theta, data, states, adj_mat=0, ll=F, marginal=T,
                       arp=1, parp=1, constraint='o', global.neighbour=T,
                       global.self=F, exo=0){
  
  # returns the complete log-likelihood, up to a constant
  # passed into optim() for the maximisation step of mzip.EM
  
  n <- length(data[, 1])
  maxar = max(arp, parp)
  
  if (constraint=='g'){
    fit <- mzip.fitted.g(theta, data, adj_mat, ll=ll, global.neighbour=global.neighbour,
                         global.self=global.self, exo=exo)
    
    # restrict coefficients to be positive to guarantee positive intensity
    if (ll==F & global.neighbour==F & any(theta[1:(n*(3+2*arp))] < 0)){
      return (-Inf)}
    if (ll==F & global.neighbour==T & any(theta[1:((2+arp)*n + arp + 1)] < 0)){
      return (-Inf)}
  } else {
    fit <- mzip.fitted.o(theta, data, ll=ll, marginal=marginal, arp=arp, parp=parp,
                         adj_mat=adj_mat)
    if (ll==F){
      z <- theta.to.params(theta, n, marginal=marginal, arp=arp, parp=parp, adj_mat=adj_mat)
      d <- z$d; A <- z$A; B <- z$B
      if (any(c(d, A, B) < 0)) return(-Inf)
    }
  }
  ints <- fit$ints; probs <- fit$probs
  
  loglikelihood <- 0
  for (t in (maxar+1):length(data[1, ])){
    if (marginal) {
      for (i in 1:n){
        if (ints[i, t] < 0){
          stop('Negative intensity')
        }
        loglikelihood <- loglikelihood + (1-states[i, t])*(data[i, t]*log(ints[i, t]+1e-60) - ints[i, t]) +
          states[i, t]*log(probs[i, t]+1e-60) + (1-states[i, t])*log(1 -probs[i, t] +1e-60)
      }
    } else{
      poiscomp <- 0
      for (i in 1:n){
        poiscomp <- poiscomp + data[i, t]*log(ints[i, t]+1e-60) - ints[i, t]
      }
      loglikelihood <- loglikelihood + (1-states[t])*poiscomp +
        states[t]*log(probs[t]+1e-60) + (1-states[t])*log(1-probs[t]+1e-60)
    }
  }
  return (loglikelihood)
}
