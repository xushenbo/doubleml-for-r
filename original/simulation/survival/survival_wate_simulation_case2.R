library(foreach)
library(doParallel)
library(mvtnorm)
library(simsurv)
library(pec)
library(DescTools)
library(zoo)
library(dplyr)
library(survival)
library(quantmod)
library(tidyr)
library(SuperLearner)
library(survSuperLearner)
library(xgboost)
# library(bartMachine)
library(KernelKnn)
library(earth)
library(fdrtool)
library(simecol)
library(MASS)
library(ggplot2)
library(data.table)
library(riskRegression)
library(mstate)

## data generating process
dgp_survival <- function(n, admin.cens, time.grid){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6)
  
  data.df$Ta0 <- simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime
  data.df$Ta1 <- simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime
  data.df$Ta <- (1-a)*data.df$Ta0 + a*data.df$Ta1
  
  data.df$C[data.df$a==0] <- simsurv(dist = "exponential", lambdas = 0.06, x = cbind(x0=1, data.df[data.df$a==0, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0.1, x1=0.4, x2=-0.7, x3=-0.4, x4=-0.5, x5=0.8, x6=-0.6))$eventtime
  data.df$C[data.df$a==1] <- simsurv(dist = "exponential", lambdas = 0.08, x = cbind(x0=1, data.df[data.df$a==1, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0, x1=0.5, x2=-0.6, x3=0.2, x4=0.6, x5=0.9, x6=-0.5))$eventtime
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Ta)/time.grid)*time.grid # , data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df,ifelse(Ta<=C, 1, 0))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  return(data.df)
}

registerDoParallel(cores=207)

nsim <- 1000
n <- 2000

time.grid <- 0.01
admin.cens <- 10
taus <- c(1,2,3,4)
cens.prob <- data.frame(matrix(NA, nrow=nsim, ncol=4))
colnames(cens.prob) <- c("tau1","tau2","tau3","tau4") # ,"tau6"
# counting.sl.lib <- c("SL.mean","SL.glm","SL.step","SL.glmnet","SL.rpartPrune","SL.gam","SL.xgboost","SL.ranger") # ,"SL.nnet","SL.kernelKnn","SL.earth" # "SL.xgboost":33s; "SL.ranger":30s; "SL.ksvm":48.19s; "SL.bartMachine":slow; "SL.polymars": not working
a.sl.lib <- c("SL.mean","SL.nnet","SL.kernelKnn","SL.rpartPrune","SL.xgboost","SL.ranger","SL.glm","SL.step","SL.gam","SL.glmnet","SL.earth") # 
event.sl.lib <- cens.sl.lib <- c("survSL.km","survSL.loglogreg","survSL.pchreg","survSL.rfsrc","survSL.coxph","survSL.expreg","survSL.weibreg") # ; "survSL.gam","survSL.pchSL"
confounders <- c("x1","x2","x3","x4","x5","x6")


# dml simulation
start.time <- Sys.time()
rmst.dml.result <- run_rmst_dml_simulation(nsim, n, admin.cens, time.grid, taus, confounders, a.sl.lib, event.sl.lib, cens.sl.lib)
end.time <- Sys.time()
end.time - start.time

save(rmst.dml.result, file="rmst_wate_dml_case2.RData")

# baseline simulation
registerDoParallel(cores=100)

nsim <- 1000
n <- 2000
nboot <- 100

start.time <- Sys.time()
rmst.baseline.result <- run_rmst_baseline_simulation(nsim, n, admin.cens, time.grid, taus, nboot)
end.time <- Sys.time()
end.time - start.time
save(rmst.baseline.result, file="rmst_wate_baseline_case2.RData")

# true value
dgp_true <- function(ntrue){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=ntrue, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(ntrue,1,0.5)
  x5 <- rbinom(ntrue,1,0.5)
  x6 <- rbinom(ntrue,1,0.5)
  ps <- 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(ntrue, 1, ps) # -2+1.2*x1+1.6*x2+1.6*x3+1.6*x4, -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  ps <- Winsorize(minval = 1e-4, maxval = 1-1e-4, ps)
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6)
  data.df$naive <- 1
  data.df$ow.tilt <- with(data.df, ps*(1-ps))
  data.df$mw.tilt <- with(data.df, pmin(ps, 1-ps))
  data.df$enw.tilt <- with(data.df, ifelse(ps==0 | ps==1, 0, -(ps*log(ps)+(1-ps)*log(1-ps))))
  
  # data.df$Tj1a0 <- simsurv(dist = "weibull", lambdas = 0.04, gammas = 0.4, x = data.frame(x0=1, x1, x2, x3, x4), betas = c(x0=0.3, x1=0.8, x2=-1.2, x3=0.2, x4=0.1))$eventtime
  # data.df$Tj1a1 <- simsurv(dist = "weibull", lambdas = 0.05, gammas = 0.5, x = data.frame(x0=1, x1, x2, x3, x4), betas = c(x0=0.2, x1=1.4, x2=-0.6, x3=0.7, x4=0.3))$eventtime
  # data.df$Tj2a0 <- simsurv(dist = "weibull", lambdas = 0.05, gammas = 0.3, x = data.frame(x0=1, x1, x2, x3, x4), betas = c(x0=0.4, x1=-1.5, x2=0.8, x3=-0.5, x4=0.2))$eventtime
  # data.df$Tj2a1 <- simsurv(dist = "weibull", lambdas = 0.04, gammas = 0.4, x = data.frame(x0=1, x1, x2, x3, x4), betas = c(x0=0.2, x1=-1.0, x2=1.5, x3=-1, x4=0.4))$eventtime
  # data.df$C <- simsurv(dist = "weibull", lambdas = 1.1, gammas = 1.5, x = data.frame(a, x0=1, x1, x2, x3, x4), betas = c(a=0.2, x0=-0.5, x1=0.4, x2=-0.7, x3=-0.5, x4=-0.5))$eventtime
  # data.df$C <- simsurv(dist = "exponential", lambdas = 0.2, x = data.frame(a, x0=1, x1, x2, x3, x4), betas = c(a=0.2, x0=0.5, x1=0.4, x2=-0.7, x3=0.5, x4=0.5))$eventtime
  
  data.df$Ta0 <- ceiling(simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime*1000)/1000
  data.df$Ta1 <- ceiling(simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime*1000)/1000
  
  data.df$Da0 <- 1
  data.df$Da1 <- 1
  
  # censor the lastest observed event
  # data.df$Da0[which.max(data.df$Ta0)] <- 0
  # data.df$Da1[which.max(data.df$Ta1)] <- 0
  
  return(data.df)
}

ntrue <- 1e+06
admin.cens <- 10
taus <- c(1,2,3,4)


start.time <- Sys.time()
large.df <- dgp_true(ntrue)
end.time <- Sys.time()
end.time - start.time

large.df$Ta0[large.df$Ta0>taus[4]+1] <- taus[4]+1
large.df$Ta1[large.df$Ta1>taus[4]+1] <- taus[4]+1


start.time <- Sys.time()
true.a0 <- data.frame(s=sort(unique(large.df$Ta0)))
true.a0$rmst.a0 <- rmst_compute_true_value(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$naive)
true.a0$rmst.ow.a0 <- rmst_compute_true_value(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$ow.tilt)
true.a0$rmst.mw.a0 <- rmst_compute_true_value(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$mw.tilt)
true.a0$rmst.enw.a0 <- rmst_compute_true_value(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$enw.tilt)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
true.a1 <- data.frame(s=sort(unique(large.df$Ta1)))
true.a1$rmst.a1 <- rmst_compute_true_value(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$naive)
true.a1$rmst.ow.a1 <- rmst_compute_true_value(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$ow.tilt)
true.a1$rmst.mw.a1 <- rmst_compute_true_value(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$mw.tilt)
true.a1$rmst.enw.a1 <- rmst_compute_true_value(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$enw.tilt)
end.time <- Sys.time()
end.time - start.time


true.value.tau1 <- rmst_find_true_value(true.a0, true.a1, taus[1])
true.value.tau2 <- rmst_find_true_value(true.a0, true.a1, taus[2])
true.value.tau3 <- rmst_find_true_value(true.a0, true.a1, taus[3])
true.value.tau4 <- rmst_find_true_value(true.a0, true.a1, taus[4])

save(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4, file="rmst_wate_true_case2.RData")

