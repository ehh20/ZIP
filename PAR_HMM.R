# Code related to the Poisson autoregression model


par.fitted <- function(theta, data, ll=F){
  
  # returns a length (length(data)) vector of fitted intensities assuming the 
  # data generating process has parameters [theta]
  
  arp <- length(theta)-2
  d <- theta[1]; a <- theta[2]; b <- theta[3:(2+arp)]
  lambda <- rep(0, length(data))
  lambda[1:arp] <- d
  if (ll){
    for (t in (arp+1):length(data)){
      lambda[t] <- d + a*lambda[t-1] + sum(b*rev(log(data[(t-arp):(t-1)]+1)))
    }
    lambda <- exp(lambda)
  } else{
    for (t in (arp+1):length(data)){
      lambda[t] <- d + a*lambda[t-1] + sum(b*rev(data[(t-arp):(t-1)]))
    }
    return(lambda)
  }
}


par.likelihood <- function(theta, data, indiv=F, ll=F){
  
  # Returns the likelihood of the data under a PAR model with parameters
  # theta = c(d, a, b_1, b_2, ...)
  
  arp <- length(theta)-2
  if (ll==F & any(theta < 0)) return (-Inf)
  likelihood <- 0
  lambdas <- par.fitted(theta, data, ll=ll)
  if (indiv){indiv_list <- rep(0, length(data)-arp)}
  for (t in (arp+1):length(data)){
    ll_t <- log(dpois(data[t], lambdas[t]))
    if (indiv) {indiv_list[t] <- ll_t}
    likelihood <- likelihood + ll_t
  }
  if (indiv) {return(indiv_list)}
  return(likelihood)
}


par.mle <- function(data, arp=1, params=rep(0.1, 2+arp), ll=F){
  # Returns the MLE for a PAR model with [arp] lagged past values. 
  # Change params to initialise at different starting values
  mle <- optim(params, fn=function(theta) -par.likelihood(theta, data, ll=ll),
               gr = function(theta) -par.gradient(theta, data))
  return (mle$par)
}


par.pr <- function(data, lambdas, arp=1){
  # return a vector of residuals based on fitted values from a PAR model
  residuals <- rep(0, length(data))
  for (t in (arp-1):length(data)){
    residuals[t] <- (data[t] - lambdas[t]) / sqrt(lambdas[t])
  }
  return (residuals[(arp+1):length(data)])
}


par.pred.log_score <- function(data, theta=0, lambdas=0, train=0, arp=1, ll=F){
  # returns the sum of log scores from a 1-step ahead forecasted density 
  if (train==0){
    train <- floor(2/3*length(data))
  }
  if (all(theta == 0) & all(lambdas==0)){
    theta <- par.mle(data[1:train], arp=arp, ll=ll)
    lambdas <- par.fitted(theta, data, ll=ll)
  } else {
    if (all(lambdas == 0)){lambdas <- par.fitted(theta, data, ll=ll)}
  }
  log_score <- 0
  for (t in (train+1):length(data)){
    log_score <- log_score - log(dpois(data[t], lambdas[t]))
  }
  return(log_score)
}


par.pred.point_diff <- function(data, theta=0, lambdas=0, train=0, arp=1, ll=ll){
  # returns the difference between the 1-step ahead point forecast and the observed values
  if (train==0){
    train <- floor(2/3*length(data))
  }
  if (all(theta == 0) & all(lambdas==0)){
    theta <- par.mle(data[1:train], arp=arp, ll=ll)
    lambdas <- par.fitted(theta, data, ll=ll)
  } else {
    if (all(lambdas == 0)) {lambdas <- par.fitted(theta, data, ll=ll)}
  }
  return ((data-lambdas)[(train+1):length(data)])
}


# ----------------------------------------------------
# Code related to the hidden markov model

hmm.sim <- function(d, a, b, p11=0.9, p00=0.9, t.max=1000, burnin=100,
                    theta=0.5, cont=TRUE, noise=TRUE){
  
  # returns a simulation of a hidden markov model process as described in section 2.1 
  # Inputs:
  # p11: probability of staying in the PAR state, given that the previous state was PAR
  # p00: probability of staying in the Poi(theta) state, given that the previous state was Poi(theta) 
  # cont=TRUE means pick up previous count from last state 1 values
  # noise=TRUE means that the zero state model has a poi(theta) dist, otherwise zero. 
  
  states <- c()
  state <- 1
  l.prev <- 0
  count.prev <- 0
  values <- c()
  for (t in 1:(t.max+burnin)){
    u <- runif(1)
    if (state==1 & u>p11){
      state <- 0
    }
    else{
      if (state==0 & u>p00){
        state <- 1
      }
    }
    l <- d + a*l.prev + b*count.prev
    poiscount <- rpois(1, l)
    if (state == 1){
      values <- c(values, poiscount)
      count.prev <- poiscount
    }
    else{
      if (noise == TRUE){
        zsval <- rpois(1, theta)
      }
      else{
        zsval <- 0
      }
      values <- c(values, zsval)
      if (cont==FALSE){
        count.prev <- zsval
      }
    }
    l.prev <- l
    states <- c(states, state)
  }
  return(values[(burnin+1):length(values)])
}



hmm.baum_welch <- function(data){
  
  # Baum Welch algorithm for fitting the PAR hidden markov model to [data]
  # returns a vector of (d, a, b, gamma, [transition probabilities by row])
  
  t.len <- length(data)
  state.dist <- c(0.5, 0.5)
  tm <- matrix(c(0.85, 0.15, 0.15, 0.85), nrow=2, byrow=T) #p00, p01, p10, p11
  initparams <- par.mle(data)
  d <- initparams[1]
  a <- initparams[2]
  b <- initparams[3]
  theta <- 0.1
  params.prev <- c(d, a, b, theta, c(tm))
  params <- rep(0, 8)
  iters <- 0
  
  while(norm(params.prev - params, type="2") > 1e-5 & iters<50){
    params.prev <- params
    iters <- iters + 1
    
    # 1 expectation step
    # 1.1 forward procedure
    alpha <- matrix(0, 2, t.len)
    alpha[1, 1] <- state.dist[1]*dpois(data[1], theta)
    alpha[2, 1] <- state.dist[2]*dpois(data[1], d)
    c <- c(1/(alpha[1, 1]+alpha[2, 1]))
    alpha[, 1] <- c[1] * alpha[, 1]
    l.prev <- d
    ycx <- matrix(0, 2, t.len)
    for (t in 2:t.len){
      l <- d + a*l.prev + b*data[t-1]
      ycx[, t] <- c(dpois(data[t], theta), dpois(data[t], l))
      l.prev <- l
      alpha[, t] <- ycx[, t] * (t(tm) %*% alpha[, t-1])
      c <- c(c, 1/(alpha[1, t]+alpha[2, t]))
      alpha[, t] <- c[t] * alpha[, t]
    }
    
    # 1.2 backward procedure
    beta <- matrix(0, 2, t.len)
    beta[, t.len] <- c[t.len] * c(1, 1)
    for (t in (t.len-1):1){
      beta[, t] <- t(t(tm)*ycx[, t+1]) %*% beta[, t+1]
      beta[, t] <- c[t]*beta[, t]
    }
    
    # 1.3 update expectations
    xcy <- alpha*beta
    xcy <-  t(t(xcy)/colSums(xcy))  # normalise
    
    tcy <- matrix(0, 4, t.len-1) # column: p00, p10, p01, p11
    for (t in 1:(t.len-1)){
      b.ep <- beta[, t+1] * ycx[, t+1]
      tcyt <- (alpha[, t] %o% b.ep) * tm
      tcy[, t] <- c(tcyt)
    }
    tcy <- t(t(tcy)/colSums(tcy))  # normalise
    
    # 2. maximisation step
    # 2.1 update initial state and transition probabilities
    state.dist <- xcy[, 1]
    tm <- matrix(rowSums(tcy), nrow=2, byrow=F) / rowSums(xcy[, 1:(t.len-1)])
    
    # 2.2 MLE for model parameters
    optparams <- optim(c(d, a, b, theta),
                       fn = function(params) hmm.baum_welch_Mstep(params, data, xcy),
                       control=list(fnscale=-1))
    likelihood <- optparams$value
    newparams <- optparams$par
    d <- newparams[1]
    a <- newparams[2]
    b <- newparams[3]
    theta <- newparams[4]
    params <- c(d, a, b, theta, c(tm))
  }
  return (params)
}


hmm.baum_welch_Mstep <- function(params, y, xcy){
  if (any(params < 0)) return (-Inf) # force parameters to be positive
  d <- params[1]
  a <- params[2]
  b <- params[3]
  theta <- params[4]
  l.prev <- d
  obj <- xcy[1, 1]*log(dpois(y[1], theta)) + xcy[2, 1]*log(dpois(y[1], l.prev))
  for (t in 2:length(y)){
    l <- d + a*l.prev + b*y[t-1]
    obj <- obj + xcy[1,t]*log(dpois(y[t],theta)) + xcy[2,t]*log(dpois(y[t],l))
    l.prev <- l
  }
  return (obj)
}