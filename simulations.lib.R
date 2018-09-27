library(actuar)
set.seed(1)
randomSample <- function(c){
  
  S <- 0
  n <- 0
  
  while (S <= c){
    
    S <- S + rpareto(1, shape = 1, scale = 2)
    n <- n + 1
  }
  
  return(n)
}
  
CMCSim2 <- function(n = 1000, K = 1000, c = 800){
  
  theta <- rep(0, K)
  for (k in 1:K){
    P <- rep(0, n)
    for (j in 1:n){
      
      m <- randomSample(c)
      P[j] <- ppois(m - 1, 5, lower.tail = FALSE)
      
    }
    theta[k] <- mean(P)
    
  }
  
  return(theta)
  
}

#' @title compCMonteCarlo
#'
#' @param shape shape parameters of the Pareto distributions
#' @param scale scale parameters of the Pareto distributions
#' @param sampRow one row containing N samples
#' @param c the right limit threshold
#' @return the CMC estimate Z
#' @examples
#' compCMonteCarlo(c(1,1.1), c(2,2.2), c(20, 25), 21)
compCMonteCarlo <- function(shape, scale, sampRow, c){
  
  z <- 0
  S <- sum(sampRow)
  m <- max(sampRow)
  
  for (i in 1:length(shape)){
    
    z <- z + ppareto(max(c - sum(sampRow[-i]), max(sampRow[-i])), shape[i], scale[i], lower.tail = F)
    
  }
  
  return(z)
}



#' @title CMC.sim
#' 
#' @param K number of Zbars
#' @param n number of samples drawn
#' @param N number of factors
#' @param c right tail threshold
#' @param parPars a dataframe with the first column the shape and the second column the 
#' scale of the Pareto distributions.
#' @return list of CMC estimate, Monte Carlo estimate and mean CMC
#' @examples
#' CMC.sim(n = 2000, N = 10, c = 100, parPars = NULL)
CMCSim <- function(n = 200, N = 10, c = 100, parPars = NULL, K = 100){
  
  
  if (is.null(parPars)){
    shapeV <- seq(.8, 1.2, length.out = N)
    scaleV <- seq(1.8, 2.2, length.out = N)
  
    parPars <- data.frame(shapeV, scaleV) 
  }
  
  Zbar <- rep(0,K)
  muBar <- rep(0,K)
  
  for (k in 1:K){
  
    dat <- data.frame(apply(parPars, 1, function(x){return(rpareto(n, x[1], x[2]))}))
  
    ZSum <- 0
    muHatSum <- 0
  
    for (i in 1:n){
    
      ZSum <- ZSum + compCMonteCarlo(shapeV, scaleV, dat[i,], c)
      muHatSum <- muHatSum + as.numeric(sum(dat[i,]) > c)
    }
  
    Zbar[k] <- ZSum/n
    muBar[k] <- muHatSum/n
    
  }
  
  return(list("Zbar" = Zbar, "muBar" = muBar))
  
}

#' @title simErrorRatio
#'
#'@description varying c  
#' @param n vector of sample size
#' @param N number of factors
#' @param c right tail threshold
#' @param parPars a dataframe with the first column the shape and the second column the 
#' scale of the Pareto distributions.
#' @param kappa the 
#' @return list of CMC estimate, Monte Carlo estimate and mean CMC
#' @examples
#' CMC.sim(n = 2000, N = 10, c = 100, parPars = NULL)
simErrorRatio <- function(n = 2^(seq(1,10)), N = 10, c = 100, K = 100, kappa = .05){
  
  shapeV <- seq(.8, 1.2, length.out = N)
  scaleV <- seq(1.8, 2.2, length.out = N)
  parPars <- data.frame(shapeV, scaleV) 
  
  print("Estimating mu")
  
  mu <- CMCSim(1e5, N, c, parPars, K = 1)$Zbar
  
  cmc <- rep(0, length(n))
  mc <- rep(0, length(n))
  PRatio <- rep(0, length(n))
  for (i in 1:length(n)){
    print(paste("n = ", n[i]))
    tmp <- CMCSim(n[i], N, c, parPars, K)
    cmc[i] <- mean(tmp$Zbar)
    mc[i] <- mean(tmp$muBar)
    
    PRatio[i] <- log(sum(abs(tmp$Zbar - mu) > kappa * mu)/K)/(sum(abs(tmp$muBar - mu) > kappa * mu)/K)
  }
  
  return(list("PRatio"= PRatio, "n" = n, "mu" = mu, "cmc" = cmc, "mc" = mc, "N" = N, "c" = c))
  
}





