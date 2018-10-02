require(ggplot2)
source("simulations.lib.R")

out <- simErrorRatio(n = seq(1, 1501, length.out = 31), kappa = .01) 

ggplot() + geom_point(aes(x = out2s$n, y = out2s$PRatio)) + geom_smooth(aes(x = out2s$n, y = out2s$PRatio), method = "lm")

ggplot() + geom_point(aes(x = out2$n, y = out2$PRatio)) + geom_smooth(aes(x = out2$n, y = out2$PRatio), method = "lm")

summary(lm(out2$PRatio ~ out2$n))


library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
summary(lm(out2s$PRatio ~ out2s$n))