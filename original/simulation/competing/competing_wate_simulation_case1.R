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

# data generating process
dgp_competing <- function(n, admin.cens, time.grid){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(0.3+0.2*x1+0.3*x2+0.3*x3-0.2*x4-0.3*x5-0.2*x6))) # 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6)
  
  data.df$Tj1a0 <- simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime
  data.df$Tj1a1 <- simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime
  data.df$Tj1 <- (1-a)*data.df$Tj1a0 + a*data.df$Tj1a1
  
  data.df$Tj2a0 <- simsurv(dist = "exponential", lambdas = 0.1, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.12, x1=-0.1, x2=0.3, x3=0.1, x4=0.2, x5=-0.4, x6=0.5))$eventtime
  data.df$Tj2a1 <- simsurv(dist = "exponential", lambdas = 0.08, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=-0.2, x2=-0.1, x3=0.2, x4=0.3, x5=0.3, x6=-0.3))$eventtime
  data.df$Tj2 <- (1-a)*data.df$Tj2a0 + a*data.df$Tj2a1
  
  data.df$C[data.df$a==0] <- simsurv(dist = "exponential", lambdas = 0.12, x = cbind(x0=1, data.df[data.df$a==0, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0.1, x1=0.4, x2=-0.7, x3=-0.4, x4=-0.5, x5=0.8, x6=-0.6))$eventtime
  data.df$C[data.df$a==1] <- simsurv(dist = "exponential", lambdas = 0.14, x = cbind(x0=1, data.df[data.df$a==1, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0, x1=0.5, x2=-0.6, x3=0.2, x4=0.6, x5=0.9, x6=-0.5))$eventtime
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Tj1, data.df$Tj2)/time.grid)*time.grid # , data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df, ifelse(pmin(Tj1,Tj2)>C,0, ifelse(Tj1<Tj2,1,2)))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  return(data.df)
}

# dml simulation
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

start.time <- Sys.time()
rmtlj.dml.result <- run_rmtlj_dml_simulation(nsim, n, admin.cens, time.grid, taus, confounders, a.sl.lib, event.sl.lib, cens.sl.lib)
end.time <- Sys.time()
end.time - start.time

save(rmtlj.dml.result, file="rmtlj_wate_dml_case1.RData")

# baseline simulation
registerDoParallel(cores=100)

nsim <- 1000
n <- 2000
nboot <- 100

time.grid <- 0.01
admin.cens <- 10
taus <- c(1,2,3,4)
cens.prob <- data.frame(matrix(NA, nrow=nsim, ncol=4))
colnames(cens.prob) <- c("tau1","tau2","tau3","tau4")
confounders <- c("x1","x2","x3","x4","x5","x6")

start.time <- Sys.time()
rmtlj.baseline.result <- run_rmtlj_baseline_simulation(nsim, n, admin.cens, time.grid, taus, nboot)
end.time <- Sys.time()
end.time - start.time

save(rmtlj.baseline.result, file="rmtlj_wate_baseline_case1.RData")

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
  ps <- 1/(1+exp(-(0.3+0.2*x1+0.3*x2+0.3*x3-0.2*x4-0.3*x5-0.2*x6))) # -1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6
  a <- rbinom(ntrue, 1, ps) # -2+1.2*x1+1.6*x2+1.6*x3+1.6*x4, -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  # ps <- Winsorize(minval = 1e-4, maxval = 1-1e-4, ps)
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6)
  data.df$naive <- 1
  data.df$ow.tilt <- with(data.df, ps*(1-ps))
  data.df$mw.tilt <- with(data.df, pmin(ps, 1-ps))
  data.df$enw.tilt <- with(data.df, ifelse(ps==0 | ps==1, 0, -(ps*log(ps)+(1-ps)*log(1-ps))))
  
  data.df$Tj1a0 <- ceiling(simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime*1000)/1000
  data.df$Tj1a1 <- ceiling(simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime*1000)/1000
  data.df$Tj2a0 <- ceiling(simsurv(dist = "exponential", lambdas = 0.1, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.12, x1=-0.1, x2=0.3, x3=0.1, x4=0.2, x5=-0.4, x6=0.5))$eventtime*1000)/1000
  data.df$Tj2a1 <- ceiling(simsurv(dist = "exponential", lambdas = 0.08, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=-0.2, x2=-0.1, x3=0.2, x4=0.3, x5=0.3, x6=-0.3))$eventtime*1000)/1000
  
  data.df$Ta0 <- pmin(data.df$Tj1a0,data.df$Tj2a0) # round(pmin(data.df$Tj1a0,data.df$Tj2a0),digits=3)
  data.df$Ta1 <- pmin(data.df$Tj1a1,data.df$Tj2a1) # round(pmin(data.df$Tj1a1,data.df$Tj2a1),digits=3)
  data.df$Da0 <- with(data.df,ifelse(Tj1a0<Tj2a0,1,2))
  data.df$Da1 <- with(data.df,ifelse(Tj1a1<Tj2a1,1,2))
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
true.a0$rmtlj.a0 <- compute_true(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$naive)
true.a0$rmtlj.ow.a0 <- compute_true(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$ow.tilt)
true.a0$rmtlj.mw.a0 <- compute_true(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$mw.tilt)
true.a0$rmtlj.enw.a0 <- compute_true(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$enw.tilt)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
true.a1 <- data.frame(s=sort(unique(large.df$Ta1)))
true.a1$rmtlj.a1 <- compute_true(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$naive)
true.a1$rmtlj.ow.a1 <- compute_true(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$ow.tilt)
true.a1$rmtlj.mw.a1 <- compute_true(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$mw.tilt)
true.a1$rmtlj.enw.a1 <- compute_true(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$enw.tilt)
end.time <- Sys.time()
end.time - start.time


true.value.tau1 <- find_value(true.a0, true.a1, taus[1])
true.value.tau2 <- find_value(true.a0, true.a1, taus[2])
true.value.tau3 <- find_value(true.a0, true.a1, taus[3])
true.value.tau4 <- find_value(true.a0, true.a1, taus[4])

save(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4, file="rmtlj_wate_true_case1.RData")

