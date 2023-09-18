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


zip.pred.density <- function(theta, data, h=1, nsims=1000, last_intensity=0,
                             ll=F, arp=1, parp=1, ret.last_intensity=F){
  
  # returns a h x (s+1) matrix of bootstrapped probabilities, where s is the maximum
  # value that is obtained in the bootstrapped. The i, j-th entry of the matrix
  # is the proportion of simulations that result in a value of j-1 on the i-th
  # 'observation' after data.
  
  sims <- matrix(0, nsims, h)
  
  if (last_intensity==0){
    fitted <- zip.fitted(theta, data, ll=ll, arp=arp, parp=parp)
    last_intensity <- fitted$ints[length(fitted$ints)]
  }
  maxar = max(arp, parp)
  past_vals = data[(length(data)-maxar+1):length(data)]
  
  for (i in 1:nsims){
    sims[i, ] <- zip.sim(theta, arp=arp, parp=parp, t.max=h, ll=ll, burnin=0,
                         past_vals=past_vals, past_intensity=last_intensity)
  }
  
  max_sim <- max(sims)
  densities_matrix <- matrix(0, max_sim+1, h)
  
  for (t in 1:h){
    h_table <- table(sims[, t])
    for (i in 0:max_sim){
      index <- match(i, names(h_table), nomatch=-1)
      if (index != -1){
        densities_matrix[i+1, t] <- as.numeric(h_table)[index]
      }
    }
  }
  
  densities_matrix <- densities_matrix / nsims
  if (ret.last_intensity){
    return(list("li"=last_intensity, "dm"=densities_matrix))
  }
  return(densities_matrix)
}


zip.pred.PIT <- function(data, ints=0, probs=0, train=0, h=1, ll=F, arp=1, parp=1,
                         ret.log_score=F, refit=F, refit_freq=0, Czado=F, J=20){
  
  # function to assess predictive performance of the univariate model. 
  # Inputs:
  # ints and probs are optional arguments. If not provided, will use the MLE obtained
  # from the first [train] timesteps. 
  # refit_freq: an integer for how often the model should be refitted. 
  # J: integer for the number of bins used for Czado's nonrandomised mean quantile function. 
  # Returns:
  # if ret.log_score=T: returns only the summed log scores
  # if Czado=T: returns a list containing the summed log scores, 
  # a length [length(data) - train] vector of quantiles,
  # and a vector of Czado's function evaluated at 1/J, 2/J, ...
  # otherwise: returns the vector of quantiles
  
  
  if (train==0){
    train <- floor(length(data)*2/3)
  }
  
  if (all(ints==0)){
    theta <- zip.EM(data[1:train], ll=ll, arp=arp, parp=parp)
    fitted <- zip.fitted(theta, data, arp=arp, parp=parp)
    ints <- fitted$ints; probs <- fitted$probs
  }
  
  PIT <- rep(0, (length(data)-train))
  if (Czado){
    lower <- rep(0, length(data)-train)
    upper <- rep(0, length(data)-train)
  }
  
  log_score <- 0
  refit_counter <- 1
  
  out_of_bounds <- c()
  
  for (t in (train+1):length(data)){
    if (refit==T & refit_counter == refit_freq){
      refit_counter <- 0
      theta <- zip.EM(data[1:(t-1)], ll=ll, arp=arp, parp=parp)
      last <- zip.fitted(theta, data[1:(t-1)], ll=ll, arp=arp, parp=parp)
      ints <- fitted$ints; probs <- fitted$probs
    }
    refit_counter <- refit_counter + 1
    
    if (h>1){
      # bootstrapped densities for more than 1 step ahead
      est_densities <- zip.pred.density(theta, data[1:(t-1)], h=h, last_intensity=ints[t],
                                        ll=ll, arp=arp, parp=parp, nsims=1e5)[, h]
      if (data[t] >= length(est_densities)){
        l = 1; u = 1
        log_score <- log_score + 10  #??
        out_of_bounds <- c(out_of_bounds, t)
      } else {
        if (data[t] == 0){
          l = 0
        } else {
          l = sum(est_densities[1:(data[t])])
        }
        u <- l + est_densities[data[t]+1]
        log_score <- log_score - log(est_densities[data[t]+1])
      }
    }
    
    if (h==1){
      # 1-step ahead forecast known for current distribution
      if (data[t] == 0){
        l = 0
        log_score <- log_score - log(probs[t] + (1-probs[t])*exp(-ints[t]))
      } else{
        l = probs[t] + (1-probs[t])*ppois(data[t]-1, ints[t])
        log_score <- log_score + -log((1-probs[t])*dpois(data[t], ints[t]))
      }
      u = probs[t] + (1-probs[t])*ppois(data[t], ints[t])
    }
    
    # randomised to get non-discrete quantiles
    uniform <- runif(1)
    if (Czado){
      lower[t-train] <- l
      upper[t-train] <- u
    }
    PIT[t-train] <- (1-uniform)*l + uniform*u
  }
  if (h>1){
    cat("Out of bootstrapped bounds:", out_of_bounds)
  }
  if (Czado){
    Czado_PIT <- c()
    for (j in 1:J){
      cdf <- 0
      for (t in 1:(length(data)-train)){
        if(j/J > lower[t]){
          if (j/J < upper[t]){
            cdf <- cdf + (j/J - lower[t]) / (upper[t] - lower[t])
          } else {
            cdf <- cdf + 1
          }
        }
      }
      Czado_PIT <- c(Czado_PIT, cdf/(length(data)-train))
    }
  }
  if (ret.log_score){
    return (log_score)
  }
  if (Czado){
    return (list('log_score'=log_score, 'PIT'=PIT, 'Czado'=Czado_PIT))
  } else{
    return (PIT)
  }
}


zip.pred.point_diff <- function(data, theta=0, train=0, h=1, ll=F, arp=1, parp=1){
  
  # returns a length [length(data) - train] vector of differences between the 
  # observed values and the h-steps ahead point forecasts. 
  # theta can optionally be specified. If not provided, will use the MLE obtained
  # from the first [train] timesteps. 
  
  if (train==0){
    train <- floor(length(data)*2/3)
  }
  
  if (all(theta==0)){
    theta <- zip.EM(data[1:train], ll=ll, arp=arp, parp=parp)
    fitted <- zip.fitted(theta, data, arp=arp, parp=parp)
    ints <- fitted$ints; probs <- fitted$probs
  }
  
  difference <- rep(0, length(data)-train-h+1)
  for (t in (train+h):length(data)){
    expected <- zip.sim(theta, arp=arp, parp=parp, t.max=h, burnin=0,
                        ll=ll, expected=T, past_intensity=ints[t-h],
                        past_vals=data[(t-arp-h+1):(t-h)])
    difference[t-train-h+1] <- data[t] - expected[h]
  }
  return(difference)
}


zip.pred.confidence_interval <- function(data, h=1, p=0.9, ll=F, arp=1, parp=1, nsims=1000){
  
  # returns a list containing a length h vector of 1, 2,..., h-step ahead forecasts, 
  # and a 2 x h matrix with the lower and upper bounds of the p*100% confidence intervals 
  # (without Bonferroni correction). CIs are as symmetric as possible and not the shortest. 
  
  theta <- zip.EM(data, ll=ll, arp=arp, parp=parp)
  maxar <- max(arp, parp)
  t.len <- length(data)
  pd <- zip.pred.density(theta, data, h=h, ll=ll, arp=arp, parp=parp,
                         ret.last_intensity=T, nsims=nsims)
  point_forecast <- zip.sim(theta, arp=arp, parp=parp, t.max=h, burnin=0, ll=F,
                            forecast=T, past_vals=data[(t.len-arp+1):t.len],
                            intensity=pd$li, expected=T)
  
  intervals <- matrix(0, 2, h)
  for (t in 1:h){
    lower <- floor(point_forecast[t])
    upper <- ceiling(point_forecast[t])
    cover <- pd$dm[lower+1, t]
    iter <- 0
    while (cover < 0.9){
      iter <- iter + 1
      if (iter %% 2 == 1){
        if (lower > 0){
          cover <- cover + pd$dm[lower, t]
          lower <- lower - 1
        }
      } else {
        if (upper < length(pd$dm[, 1])){
          upper <- upper + 1
          cover <- cover + pd$dm[upper, t]
        }
      }
    }
    intervals[1, t] <- lower
    intervals[2, t] <- upper
  }
  return (list("forecast"=point_forecast, "CI"=intervals))
}


zip.pr <- function(data, ints, probs, maxar=1){
  # return a vector of pearson residuals based on fitted intensities [ints] and 
  # probabilities [probs] without those from the first [maxar] observations
  
  residuals <- rep(0, length(data)-maxar)
  for (t in (maxar+1):(length(data))){
    m <- (1-probs[t]) * ints[t]
    v <- m * (1 + probs[t]*ints[t])
    residuals[t] <- (data[t] - m) / sqrt(v)
  }
  return (residuals)
}


zip.rqr <- function(data, ints, probs){
  # returns randomised quantile residuals
  
  residuals <- rep(0, (length(data)-1))
  for (t in 1:(length(data)-1)){
    if (data[t+1] == 0){
      a = 0
    } else{
      a = probs[t] + (1-probs[t])*ppois(data[t+1]-1, ints[t])
    }
    b = probs[t] + (1-probs[t])*ppois(data[t+1], ints[t])
    residuals[t] <- runif(1, min=a, max=b)
  }
  return (residuals)
}


vuong_nonnested <- function(lls1, lls2){
  
  # inputs:
  # lls1: length T vector of log likelihoods under the first model
  # lls2: length T vector of log likelihoods under the second model
  # returns the test statistic for Vuong's test for non-nested models
  
  n <- length(lls1)
  E <- 0
  w <- 0
  for (t in 1:n){
    E <- E + lls1[t] - lls2[t]
    w <- w + (lls1[t] - lls2[t])^2
  }
  E <- E / sqrt(n)
  w <- w/n - E^2
  return(E/w)
}
