library(actuar)
set.seed(1)
random.sample <- function(c){
  
  S <- 0
  n <- 0
  
  while (S <= c){
    
    S <- S + rpareto(1, shape = 1, scale = 2)
    n <- n + 1
  }
  
  return(n)
}

CMC.sim2 <- function(n = 1000, K = 1000, c = 800){
  
  theta <- rep(0, K)
  for (k in 1:K){
    P <- rep(0, n)
    for (j in 1:n){
      
      m <- random.sample(c)
      P[j] <- ppois(m - 1, 5, lower.tail = FALSE)
      
    }
    theta[k] <- mean(P)
    
  }
  
  return(theta)
  
}

#' @title CMC.sim
#'
#' @param n number of samples drawn
#' @param N number of factors
#' @param c right tail threshold
#' @param parPars a dataframe with the first column the shape and the second column the 
#' scale of the Pareto distributions.
#' @return list of CMC estimate, Monte Carlo estimate and mean CMC
#' @examples
#' CMC.sim(n = 2000, N = 10, c = 100, parPars = NULL)
CMC.sim <- function(n = 2000, N = 10, c = 100, parPars = NULL){
  
  
  if (is.null(parPars)){
    shapeV <- seq(.8, 1.2, length.out = N)
    scaleV <- seq(1.8, 2.2, length.out = N)
  
    parPars <- data.frame(shapeV, scaleV) 
  }
  
  dat <- data.frame(apply(parPars, 1, function(x){return(rpareto(n, x[1], x[2]))}))
  
  Z <- rep(0, n)
  muHat <- 0
  
  for (i in 1:n){
    
    Z[i] <- compCMonteCarlo(shapeV, scaleV, dat[i,], c)
    muHat <- muHat + as.numeric(sum(dat[i,]) > c)
  }
  
  return(list("Z" = Z, "muHat" = muHat/n, "Zbar" = mean(Z)))
  
}

compCMonteCarlo <- function(shape, scale, sampRow, c){
  
  z <- 0
  S <- sum(sampRow)
  m <- max(sampRow)
  
  for (i in 1:length(shape)){
    
    z <- z + ppareto(max(c - sum(sampRow[-i]), max(sampRow[-i])), shape[i], scale[i], lower.tail = F)
    
  }
  
  return(z)
}

