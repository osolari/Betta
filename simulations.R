require(ggplot2)
source("simulations.lib.R")

out2s <- simErrorRatio(kappa = .01) 

ggplot() + geom_point(aes(x = out$n, y = out$PRatio)) + geom_smooth(aes(x = out$n, y = out$PRatio), method = "lm")


ggplot() + geom_point(aes(x = out2$n, y = out2$PRatio)) + geom_smooth(aes(x = out2$n, y = out2$PRatio), method = "lm")

summary(lm(out2$PRatio ~ out2$n))


library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
