zip.sim <- function(theta, arp=1, parp=1, ll=F, t.max=500, burnin=100,
                    past_vals=0, past_intensity=0, expected=F,
                    zero_inflated=T){
  
  # returns a length t.max ZIP process
  # set zero_inflated=F to return a PAR(1, q) process
  # set expected=T and input past_vals and past_intensity to return a length
  # t.max-steps ahead point forecast
  
  d <- theta[1]; a <- theta[2]; b <- theta[3:(2+arp)]
  if (zero_inflated){
    gamma <- theta[3+arp]; alpha <- theta[4+arp]
    beta <- theta[(5 + arp):(4 + arp + parp)]
  }
  maxar = max(arp, parp)
  
  if (burnin==0 & length(past_vals)>=maxar & past_intensity > 0){
    # for forecasting
    values <- c(past_vals, rep(0, t.max))
    l.prev <- past_intensity; mu.prev <- log(past_intensity)
  } else{
    values <- rep(0, t.max + burnin + maxar)
    l.prev <- d; mu.prev <- d
  }
  
  for (t in (maxar + 1):(t.max + burnin + maxar)){
    if (ll==T){
      mu <- d + a*mu.prev + sum(b*rev(log(values[(t-arp):(t-1)]+1)))
      l <- exp(mu)
      mu.prev <- mu
    } else{
      l <- d + a*l.prev + sum(b*rev(values[(t-arp):(t-1)]))
    }
    if (zero_inflated){
      p <- 1 / (1 + exp(-(gamma + alpha*l.prev + sum(beta*rev(values[(t-parp):(t-1)])))))
    } else {
      p <- 0
    }
    if (expected){
      values[t] <- l*(1-p)
    } else{
      if (runif(1) < p){
        values[t] <- 0
      } else{
        values[t] <- rpois(1, l)}
    }
    l.prev <- l
  }
  return(values[(burnin+maxar+1):length(values)])
}


zip.fitted <- function(theta, data, ll=F, arp=1, parp=1, ret.forecast=F, exo=0){
  
  # returns a list containing a vector of fitted intensities and a vector of fitted
  # zero-inflation probabilities
  
  # read in parameters
  d <- theta[1]; a <- theta[2]; b <- theta[3:(2+arp)]
  gamma <- theta[3+arp]; alpha <- theta[4+arp]
  beta <- theta[(5 + arp):(4 + arp + parp)]
  if (any(exo!=0)){
    c <- theta[5 + arp + parp]
  }
  maxar <- max(arp, parp)
  
  if (ret.forecast){
    end <- length(data) + 1
  } else{
    end <- length(data)
  }
  intensities <- rep(0, end)
  probs <- rep(0, end)
  intensities[maxar] <- theta[1]
  
  # recursions for obtaining next fitted values
  for (t in (maxar+1):end){
    if (ll==T){
      if (any(exo!=0)){
        intensities[t] <- d + a*intensities[t-1] + sum(b*rev(log(data[(t-arp):(t-1)]+1))) + c*exo[t]
      } else {
        intensities[t] <- d + a*intensities[t-1] + sum(b*rev(log(data[(t-arp):(t-1)]+1)))
      }
      probs[t] <- 1 / (1 + exp(-(gamma + alpha*exp(intensities[t-1]) + sum(beta*rev(data[(t-parp):(t-1)])))))
    }
    else {
      intensities[t] <- d + a*intensities[t-1] + sum(b*rev(data[(t-arp):(t-1)]))
      probs[t] <- 1 / (1 + exp(-(gamma + alpha*intensities[t-1] + sum(beta*rev(data[(t-parp):(t-1)])))))
    }
  }
  if (ll==T){
    intensities <- exp(intensities)
  }
  if(ret.forecast){
    return (list(int=intensities[end], prob=probs[end]))
  }
  return (list(ints=intensities, probs=probs))
}


zip.PL <- function(theta, data, ll=F, indiv=F, arp=1, parp=1, exo=0){
  
  # returns the (incomplete) log-likelihood for [data] for ZIP model with parameters [theta]
  
  loglikelihood <- 0
  maxar <- max(arp, parp)
  if (indiv) {indiv.l <- rep(0, length(data)-maxar)}
  f <- zip.fitted(theta, data, ll=ll, arp=arp, parp=parp, exo=exo)
  
  for (t in (maxar + 1):length(data)){
    llt <- log(f$probs[t]*(data[t]==0) + (1-f$probs[t])*dpois(data[t], f$ints[t]))
    if (indiv) {indiv.l[t-1] <- llt}
    loglikelihood <- loglikelihood + llt
  }
  if (indiv) {return(indiv.l)}
  return (loglikelihood)
}


zip.EM <- function(data, ll=F, arp=1, parp=1, lchange=1e-5, max.iters=50, exo=0){
  
  # Expectation-Maximisation algorithm for fitting a ZIP model to [data]
  
  # initialise values
  theta <- c(rep(0.1, 2+arp), rep(0, 2+parp+any(exo!=0)))
  t.len = length(data)
  iters <- 0
  maxar = max(arp, parp)
  pl.prev <- 1
  pl <- 2
  pl_list <- c()
  
  while(abs((pl.prev - pl)/pl.prev) > lchange & iters<max.iters){
    iters <- iters + 1
    pl.prev <- pl
    states <- rep(0, t.len)
    
    # Expectation step
    f <- zip.fitted(theta, data, ll=ll, arp=arp, parp=parp, exo=exo)
    for (t in (maxar + 1):t.len){
      if (data[t]==0) {
        states[t] <-  f$probs[t] / (f$probs[t] + (1-f$probs[t])*exp(-f$ints[t]))
      }else {
        states[t] <-  0}
    }
    
    #Maximisation step
    theta <- optim(theta, fn = function(theta) zip.Mstep(theta, data, states,
                                                         ll=ll, arp=arp, parp=parp, exo=exo),
                   control=list(fnscale=-1))$par
    
    pl <- zip.PL(theta, data, ll=ll, arp=arp, parp=parp, exo=exo)
    pl_list <- c(pl_list, pl)
  }
  if (iters > max.iters){
    cat("\nFailed to convege in", max.iters, "iterations")
  }
  return (theta)
}


zip.Mstep <- function(theta, data, states, ll=F, arp=1, parp=1, exo=0){
  
  # complete data log-likelihood to pass into optim() for the maximisation step
  # of the zip.EM function
  
  if(ll==F){
    if (any(theta[1:(2+arp)] <= 0)) return (-Inf)
  }
  maxar <- max(arp, parp)
  logl <- 0
  
  f <- zip.fitted(theta, data, ll=ll, arp=arp, parp=parp, exo=exo)
  l <- f$ints; p <- f$probs
  for (t in (maxar + 1):length(data)){
    # + small constant to prevent underflow in log
    logl <- logl + (1-states[t])*(data[t]*log(l[t]+1e-30) - l[t]) +
      states[t]*log(p[t]+1e-30) + (1-states[t])*log(1-p[t]+1e-30)
  }
  return (logl)
}

