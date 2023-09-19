source("ZIP.R")
source("MZIP_helper_functions.R")


mzip.sim <- function(theta, n=2, t.max=500, burnin=100, ll=F, marginal=T, arp=1,
                     parp=1, ret.fitted=F, adj_mat=0, exo=0, exo_coefs=0,
                     # past_vals=0, past_intensity=0, expected=F TODO: k>1 step ahead forecasting
){
  
  # returns a n x t.max matrix of a simulation of a multivariate ZIP process
  
  z <- theta.to.params(theta, n, adj_mat=adj_mat, marginal=marginal,
                       arp=arp, parp=parp)
  # A, B are n x n matrices
  # pd is a scalar, pa, pb are length n vectors
  # if marginal=T, pd, pa, pb have the same dimensions as d, A, B
  d <- z$d; A <- z$A; B <- z$B; pd <- z$pd; pa <- z$pa; pb <- z$pb
  
  l.prev <- rep(0, n)
  mu.prev <- rep(0, n)
  count.prev <- rep(0, n*arp)
  p.count.prev <- rep(0, n*parp)
  values <- matrix(0, n, (t.max+burnin))
  intlist <- matrix(0, n, (t.max+burnin))
  ptlist <- rep(0, n, (t.max+burnin))
  if (marginal){
    ptlist <- matrix(0, n, (t.max+burnin))
  }
  
  for (t in 1:(t.max+burnin)){
    if (ll==F){
      l <- d + A %*% l.prev + B %*% count.prev
    } else {
      if (any(exo!=0) & t > burnin){
        mu <- d + A %*% mu.prev + B %*% log(count.prev + 1) + exo_coefs*exo[, t-burnin]
      } else {
        mu <- d + A %*% mu.prev + B %*% log(count.prev + 1)
      }
      mu.prev <- mu
      l <- exp(mu)
    }
    count <- rep(0, n)
    intlist[, t] <- l
    
    if (marginal) {
      pt <- 1 / (1 + exp(-(pd + pa %*% l.prev + pb %*% p.count.prev)))
      ptlist[, t] <- pt
      for (i in 1:n){
        if (runif(1) > pt[i]) {count[i] <- rpois(1, l[i])}
      }
    } else{
      pt <- 1 / (1 + exp(-(pd + sum(pa*l.prev) + sum(pb*p.count.prev))))
      ptlist[t] <- pt
      if (runif(1) > pt){
        for (i in 1:n) {
          count[i] <- rpois(1, l[i])
        }
      }
    }
    values[, t] <- count
    l.prev <- l
    if (arp==1) {
      count.prev <- count
    } else{
      count.prev <- c(count, count.prev[1:((arp-1)*n)])
    }
    if (parp==1) {
      p.count.prev <- count
    } else{
      p.count.prev <- c(count, p.count.prev[1:((parp-1)*n)])
    }
  }
  if (ret.fitted){
    return (list('values'=values[,(burnin+1):(t.max+burnin)],
                 'ints'=intlist[,(burnin+1):(t.max+burnin)],
                 'probs'=ptlist[(burnin+1):(t.max+burnin)]))
  }
  return(values[,(burnin+1):(t.max+burnin)])
}


mzip.PL <- function(theta, data, ll=F, adj_mat=0, marginal=T,
                    arp=1, parp=1, constraint='o', exo=0, global.neighbour=T,
                    global.self=F){
  # returns the (incomplete data) log-likelihood of [data] under a multivariate 
  # ZIP model with parameters [theta]
  
  f <- mzip.fitted(theta, data, constraint=constraint, ll=ll, marginal=marginal,
                   arp=arp, parp=parp, adj_mat=adj_mat, exo=exo,
                   global.neighbour=global.neighbour, global.self=global.self)
  maxar <- max(arp, parp)
  n <- length(data[, 1])
  loglikelihood <- 0
  if(marginal){
    for (t in (1+maxar):length(data[1,])){
      for (i in 1:n){
        loglikelihood <- loglikelihood + log(f$probs[i,t]*(data[i, t]==0) +
                                               (1-f$probs[i,t])*dpois(data[i, t], f$ints[i, t]))
      }
    }
  } else {
    for (t in (1+maxar):length(data[1,])){
      poiscomp <- 1
      for (i in 1:n){
        poiscomp <- poiscomp*dpois(data[i, t], f$ints[i, t])}
      loglikelihood <- loglikelihood + log(f$probs[t]*all(data[, t]==0) + (1-f$probs[t])*poiscomp)
    }
  }
  return(loglikelihood)
}


mzip.fitted <- function(theta, data, constraint='o', ll=F, marginal=F, arp=1, parp=1, adj_mat=0,
                        exo=0, global.neighbour=T, global.self=F){
  # returns a list containing a matrix of intensities, and a matrix / vector (if marginally zero-inflated)
  # of zero-inflation probabilities, depending on the model specified in the input. 
  # note: mzip.fitted.o and mzip.fitted.g can be found in MZIP.helper_functions.R
  
  if (constraint=='g'){
    return (mzip.fitted.g(theta, data, adj_mat, ll=ll, exo=exo, global.self=global.self,
                          global.neighbour=global.neighbour))
  }
  if (constraint=='o'){
    return (mzip.fitted.o(theta, data, ll=ll, marginal=marginal, arp=arp, parp=parp, adj_mat=adj_mat))
  }
}


mzip.EM <- function(data, adj_mat=0, ll=F, arp=1, parp=1, marginal=T, max.iters=50,
                    ret.llist=F, lchange=0.0001, univar=T, check.convergence=F,
                    M_maxits=500, constraint='o', exo=0, global.neighbour=T, global.self=F){
  
  # fits a multivariate ZIP model using the Expectation-Maximisation algorithm
  # returns the list of MLEs as a concatenated vector in the order d, a, b, gamma, alpha, beta
  # optionally returns a list of log-likelihoods by setting ret.llist=T
  
  # set initial parameters
  n <- length(data[, 1])
  t.len = length(data[1, ])
  maxar <- max(arp, parp)
  iters <- 0
  llist <- c()
  pllist <- c()
  pl.prev <- 2
  pl <- 1
  converged = T
  if (all(adj_mat==0)){
    adj <- F
  } else {
    adj <- T
  }
  
  if (constraint=='g'){
    if (univar==F) {
      theta=rep(0.01, n*(6 + 2*(arp+parp)))
    } else{
      # initialised parameters by fitting univariate models
      u <- matrix(0, n, 4+arp+parp+any(exo!=0))
      for (i in 1:n){
        u[i, ] <- zip.EM(data[i, ], ll=ll, arp=arp, parp=parp, exo=exo)
      }
      if (global.neighbour){
        theta <- c(u[, 1], u[, 2], 0, u[, 3:(2+arp)], rep(0, arp),
                   u[, 3+arp], u[, 4+arp], 0, u[, (5+arp):(4+arp+parp)], rep(0, parp))
      }
      else{
        theta <- c(u[, 1], u[, 2], rep(0, n), u[, 3:(2+arp)], rep(0, n*arp),
                   u[, 3+arp], u[, 4+arp], rep(0, n), u[, (5+arp):(4+arp+parp)], rep(0, n*parp))
      }
      if (any(exo!=0)){
        theta <- c(theta, u[, 5+arp+parp])
      }
    }
  }
  
  else {
    nmat <- sum(adj_mat)
    if (adj==T & marginal==T) {lparams <- n + (1+arp)*nmat; pparams <- n + (1+parp)*nmat}
    if (adj==T & marginal==F) {lparams <- n + (1+arp)*nmat; pparams <- 2*n + 1}
    if (adj == F & marginal==T){lparams <- n + (1+arp)*n^2; pparams <- n + (1+parp)*n^2}
    if (adj==F & marginal==F) {lparams <- n + (arp+1)*n^2; pparams <- (1+parp)*n + 1}
    
    if (univar==F) {
      # initialise all parameters at 0.01
      theta=rep(0.01, (lparams+pparams))
    } else {
      # initialise by fitting univariate models
      if (adj==F){adj_mat <- matrix(1, n, n)}
      univarparams <- matrix(0, n, 4+arp+parp)
      d <- univarparams[, 1]
      A <- fill_diag(univarparams[, 2], adj_mat)
      B <- c()
      for (k in 1:arp){
        B <- c(B, fill_diag(univarparams[, 2+k], adj_mat))}
      
      if (marginal){
        pd <- univarparams[, 3+arp]
        pa <- fill_diag(univarparams[, 4+arp], adj_mat)
        pb <- c()
        for (k in 1:parp){
          pb <- c(pb, fill_diag(univarparams[, 4+arp+k], adj_mat))}
        theta <- c(d, A, B, pd, pa, pb)
      } else{
        theta <- c(d, A, B, rep(0, pparams))
      }
    }
  }
  
  while(abs((pl - pl.prev)/pl.prev) > lchange & iters<max.iters){
    
    pl.prev <- pl
    iters <- iters + 1
    fit <- mzip.fitted(theta, data, constraint=constraint, ll=ll, marginal=marginal,
                       arp=arp, parp=parp, adj_mat=adj_mat, exo=exo,
                       global.neighbour=global.neighbour, global.self=global.self)
    ints <- fit$ints; probs <- fit$probs
    
    if (marginal){
      states <- matrix(0, n, t.len)
    } else {
      states <- rep(0, t.len)
    }
    
    # Expectation step
    for (t in (maxar+1):(t.len)){
      if (marginal){
        for (i in 1:n){
          if (data[i, t] == 0){
            states[i, t] <- probs[i, t] / (probs[i, t] + (1 - probs[i, t])*exp(-ints[i, t]))
          } else {states[i, t] <- 0}
        }
      } else {
        if (all(data[, t]==0)) {
          poiscomp <- 1
          for (i in 1:n){
            poiscomp <- poiscomp*exp(-ints[i, t])}
          states[t] <-  probs[t] / (probs[t] + (1-probs[t])*poiscomp)}
        else {states[t] <-  0}
      }
    }
    
    #Maximisation step
    opt <- optim(theta, fn = function(t) mzip.Mstep(t, data, states,
                                                    adj_mat=adj_mat, ll=ll, marginal=marginal, arp=arp, parp=parp,
                                                    constraint=constraint, global.neighbour=global.neighbour,
                                                    global.self=global.self, exo=exo),
                 control=list(fnscale=-1, maxit=M_maxits))
    print(opt$counts)
    theta <- opt$par
    llist <- c(llist, opt$value)
    pl <- mzip.PL(theta, data, ll=ll, adj_mat=adj_mat,
                  marginal=marginal, arp=arp, parp=parp, constraint=constraint,
                  exo=exo, global.neighbour=global.neighbour, global.self=global.self)
    pllist <- c(pllist, pl)
  }
  if (iters==max.iters){
    print("Failed to converge")
    converged = F}
  if (ret.llist){
    return (list("theta"=theta, "llist"=llist, "pllist"=pllist))
  } else{
    if (check.convergence == T){
      return (c(theta, converged))
    } else{
      return(theta)
    }
  }
}


mzip.PIT <- function(data, fitted, train=365){
  
  # returns a length n*(T-train) vector of quantiles for hte
  # obtain the fitted values from mzip.fitted([theta],...), where [theta] is from mzip.EM on data[, 1:train]
  
  n <- length(data[, 1])
  PITs <- c()
  for (i in 1:n){
    ints <- fitted$ints[i,]; probs <- fitted$probs[i, ]
    p <- zip.pred.PIT(data[i,], ints=ints, probs=probs, h=1, ll=F, arp=1, parp=1, train=train)
    PITs <- c(PITs, p)
  }
  return(PITs)
}


mzip.Czado <- function(data, fitted, train=365){
  
  # returns a length 20 vector corresponding to Czado's mean PIT function for
  # the combined data evaluated at 0.05, 0.1, ..., 1
  
  n <- length(data[, 1])
  data.len <- length(data[1, ])
  joint_data <- c(0)
  ints <- c(0)
  probs <- c(0)
  for (i in 1:n){
    joint_data <- c(joint_data, data[i, (train+1):data.len])
    ints <- c(ints, fitted$ints[i, (train+1):data.len])
    probs <- c(probs, fitted$probs[i, (train+1):data.len])
  }
  p <- zip.pred.PIT(data[i,], ints=ints, probs=probs, h=1, ll=F, arp=1, parp=1, train=1, Czado=T)
  return(p$Czado)
}


mzip.log_score <- function(data, fitted, train=365){
  
  # returns a length n vector of the sums of log-scores for each series from the
  # h-step ahead forecasted densities
  
  n <- length(data[, 1])
  ls <- c()
  for (i in 1:n){
    ints <- fitted$ints[i,]; probs <- fitted$probs[i, ]
    p <- zip.pred.PIT(data[i,], ints=ints, probs=probs, h=1, ll=F, arp=1, parp=1, train=train, ret.log_score=T)
    ls <- c(ls, p)
  }
  return(ls)
}


mzip.point_diff_mse <- function(data, fitted, train=365){
  
  # returns a length n vector, 1 entry per series, of the mean square errors
  # between the observed values and the 1-step ahead forecast
  
  n <- length(data[, 1])
  m <- c()
  for (i in 1:n){
    d <- c()
    for (t in (train+1):length(data[1,])){
      d <- c(d, data[i, t] - (1-fitted$probs[i, t])*fitted$ints[i, t])
    }
    m <- c(m, mean(d**2))
  }
  return(m)
}