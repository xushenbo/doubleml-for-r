#
lognormal <- function(data.df, betas, sigma, n)
{return(exp(sigma*qnorm(runif(n))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

loglogistic <- function(data.df, betas, gamma, n)
{return(exp(gamma*(log(runif(n))-log(1-runif(n)))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

aft <- function(data.df, betas, sigma, n)
{return(exp(sigma*rnorm(n, mean=0, sd=1)+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

# dml
rmtlj_extract_avarage <- function(train.ueif.a1, train.ueif.a0, estimation.ueif.a1, estimation.ueif.a0, train.time, estimation.time, tau){
  ## point estimate
  # a=0
  ueif.a0 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"), all=TRUE, by="t")
  ueif.a0 <- approxTime(ueif.a0, xout=ueif.a0$t)
  if(sum(is.na(ueif.a0))>0)
  {ueif.a0 <- approxTime(rbind(0, ueif.a0), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a0 <- t(ueif.a0[,-1])
  
  estimate.a0.lcm <- gcmlcm(x=sort(unique(c(train.time, estimation.time)))[1:max(which(!is.na(colMeans(ueif.a0))))], y=colMeans(ueif.a0)[1:max(which(!is.na(colMeans(ueif.a0))))], type="gcm")[c("x.knots","y.knots")]
  names(estimate.a0.lcm) <- c("t", "estimate.a0")
  ueif.df <- merge(data.frame(t=sort(unique(c(train.time, estimation.time)))), estimate.a0.lcm, all=TRUE, by="t")
  
  # a=1
  ueif.a1 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t") , all=TRUE, by="t")
  ueif.a1 <- approxTime(ueif.a1, xout=ueif.a1$t)
  if(sum(is.na(ueif.a1))>0)
  {ueif.a1 <- approxTime(rbind(0, ueif.a1), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a1 <- t(ueif.a1[,-1])
  estimate.a1.lcm <- gcmlcm(x=sort(unique(c(train.time, estimation.time)))[1:max(which(!is.na(colMeans(ueif.a1))))], y=colMeans(ueif.a1)[1:max(which(!is.na(colMeans(ueif.a1))))], type="gcm")[c("x.knots","y.knots")]
  names(estimate.a1.lcm) <- c("t", "estimate.a1")
  ueif.df <- merge(ueif.df, estimate.a1.lcm, all=TRUE, by="t")
  
  ueif.df <- approxTime(ueif.df, xout=ueif.df$t)
  ueif.df$estimate.diff <- ueif.df$estimate.a1-ueif.df$estimate.a0
  
  ## standard error
  ueif.df$estimate.a0.se <- cummax(sqrt(colMeans((ueif.a0-matrix(ueif.df$estimate.a0, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(nrow(ueif.a0)))
  ueif.df$estimate.a1.se <- cummax(sqrt(colMeans((ueif.a1-matrix(ueif.df$estimate.a1, nrow=nrow(ueif.a1), ncol=ncol(ueif.a1), byrow=TRUE))^2))/sqrt(nrow(ueif.a1)))
  ueif.df$estimate.diff.se <- cummax(sqrt(colMeans(((ueif.a1-ueif.a0)-matrix(ueif.df$estimate.diff, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(nrow(ueif.a0)))
  
  tau.ind <- max(which((tau-ueif.df$t)>0))
  ueif.df <- ueif.df[tau.ind,]+(tau-ueif.df$t[tau.ind])*(ueif.df[tau.ind+1,]-ueif.df[tau.ind,])/(ueif.df$t[tau.ind+1]-ueif.df$t[tau.ind])
  return(ueif.df[,-1])
}

run_rmtlj_dml_simulation <- function(nsim, n, admin.cens, time.grid, taus, confounders, a.sl.lib, event.sl.lib, cens.sl.lib) {
  foreach(i=1:nsim, .combine=function(x,y) mapply(cbind, x, y, SIMPLIFY=FALSE)) %dopar% {
    sim.df <- dgp_competing(2*n, admin.cens, time.grid)
    # cens.prob[i,"tau2"] <- sum(sim.df$event==0 & sim.df$end.time < taus[1])/(2*n)
    # cens.prob[i,"tau4"] <- sum(sim.df$event==0 & sim.df$end.time < taus[2])/(2*n)
    # cens.prob[i,"tau6"] <- sum(sim.df$event==0 & sim.df$end.time < taus[3])/(2*n)
    
    train.ind <- sample(1:(2*n), n, replace=FALSE) # sample training set
    train.df <- sim.df[train.ind,]
    train.df <- train.df[order(train.df$id),]
    estimation.df <- sim.df[setdiff(1:(2*n), train.ind),]
    estimation.df <- estimation.df[order(estimation.df$id),] # estimation.df <- estimation.df[order(estimation.df$end.time,-estimation.df$event),]
    
    while(max(train.df$end.time[train.df$a==0 & train.df$event==0])<taus[4] | max(train.df$end.time[train.df$a==1 & train.df$event==0])<taus[4] |
          max(train.df$end.time[train.df$a==0 & train.df$event==1])<taus[4] | max(train.df$end.time[train.df$a==1 & train.df$event==1])<taus[4] |
          max(train.df$end.time[train.df$a==0 & train.df$event==2])<taus[4] | max(train.df$end.time[train.df$a==1 & train.df$event==2])<taus[4] |
          max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==0])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==0])<taus[4] |
          max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==1])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==1])<taus[4] |
          max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==2])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==2])<taus[4]) {
      sim.df <- dgp4_competing(2*n, admin.cens, time.grid)
      train.ind <- sample(1:(2*n), n, replace=FALSE) # sample training set
      train.df <- sim.df[train.ind,]
      train.df <- train.df[order(train.df$id),]
      estimation.df <- sim.df[setdiff(1:(2*n), train.ind),]
      estimation.df <- estimation.df[order(estimation.df$id),]}
    
    ## training set
    train.df$naive <- 1
    train.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = estimation.df$a, X = estimation.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                               cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
    train.df$iptw <- with(train.df,1/(ps*a+(1-ps)*(1-a)))
    train.df$iptw.trunc <- with(train.df,ifelse(iptw>quantile(iptw,0.99),quantile(iptw,0.99),iptw))
    train.df$ow <- with(train.df,ifelse(a,1-ps,ps))
    train.df$mw <- with(train.df,pmin(ps,1-ps)*ifelse(a,1/ps,1/(1-ps)))
    train.df$enw <- with(train.df,-(ps*log(ps)+(1-ps)*log(1-ps))*ifelse(a,1/ps,1/(1-ps)))
    train.df$ow.tilt <- with(train.df, ps*(1-ps))
    train.df$mw.tilt <- with(train.df, pmin(ps, 1-ps))
    train.df$enw.tilt <- with(train.df,-(ps*log(ps)+(1-ps)*log(1-ps)))
    
    ## estimation set
    estimation.df$naive <- 1
    estimation.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = train.df$a, X = train.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                                    cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred) # estimation.df$ps <- predict(glm(a~x1+x2+x3+x4,data=estimation.df,family=binomial()),type="response")
    estimation.df$iptw <- with(estimation.df,1/(ps*a+(1-ps)*(1-a)))
    estimation.df$iptw.trunc <- with(estimation.df,ifelse(iptw>quantile(iptw,0.99),quantile(iptw,0.99),iptw))
    estimation.df$ow <- with(estimation.df,ifelse(a,1-ps,ps))
    estimation.df$mw <- with(estimation.df,pmin(ps,1-ps)*ifelse(a,1/ps,1/(1-ps)))
    estimation.df$enw <- with(estimation.df,-(ps*log(ps)+(1-ps)*log(1-ps))*ifelse(a,1/ps,1/(1-ps)))
    estimation.df$ow.tilt <- with(estimation.df, ps*(1-ps))
    estimation.df$mw.tilt <- with(estimation.df, pmin(ps, 1-ps))
    estimation.df$enw.tilt <- with(estimation.df,-(ps*log(ps)+(1-ps)*log(1-ps)))
    
    ## train
    train.all.survival.super.a0 <- try(survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]!=0), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                        newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    while(inherits(train.all.survival.super.a0, "try-error"))
    {
      train.all.survival.super.a0 <- try(survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]!=0), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                          newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    }
    train.G.a0 <- train.all.survival.super.a0$cens.SL.predict
    train.S.a0 <- train.all.survival.super.a0$event.SL.predict
    
    train.all.survival.super.a1 <- try(survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]!=0), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                        newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    while(inherits(train.all.survival.super.a1, "try-error"))
    {
      train.all.survival.super.a1 <- try(survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]!=0), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                          newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    }
    train.G.a1 <- train.all.survival.super.a1$cens.SL.predict
    train.S.a1 <- train.all.survival.super.a1$event.SL.predict
    
    train.Sj.super.a0 <- try(survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]==1), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                              newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    while(inherits(train.Sj.super.a0, "try-error"))
    {
      train.Sj.super.a0 <- try(survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]==1), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    }
    train.Sj.a0 <- train.Sj.super.a0$event.SL.predict
    
    train.Sj.super.a1 <- try(survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]==1), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                              newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    while(inherits(train.Sj.super.a1, "try-error"))
    {
      train.Sj.super.a1 <- try(survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]==1), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    }
    train.Sj.a1 <- train.Sj.super.a1$event.SL.predict
    
    ## estimation
    estimation.all.survival.super.a0 <- try(survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]!=0), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                             newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    while(inherits(estimation.all.survival.super.a0, "try-error"))
    {
      estimation.all.survival.super.a0 <- try(survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]!=0), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                               newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    }
    estimation.G.a0 <- estimation.all.survival.super.a0$cens.SL.predict
    estimation.S.a0 <- estimation.all.survival.super.a0$event.SL.predict
    
    estimation.all.survival.super.a1 <- try(survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]!=0), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                             newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    while(inherits(estimation.all.survival.super.a1, "try-error"))
    {
      estimation.all.survival.super.a1 <- try(survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]!=0), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                               newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    }
    estimation.G.a1 <- estimation.all.survival.super.a1$cens.SL.predict
    estimation.S.a1 <- estimation.all.survival.super.a1$event.SL.predict
    
    estimation.Sj.super.a0 <- try(survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]==1), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                   newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    while(inherits(estimation.Sj.super.a0, "try-error"))
    {
      estimation.Sj.super.a0 <- try(survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]==1), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                     newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    }
    estimation.Sj.a0 <- estimation.Sj.super.a0$event.SL.predict
    
    estimation.Sj.super.a1 <- try(survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]==1), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                   newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    while(inherits(estimation.Sj.super.a1, "try-error"))
    {
      estimation.Sj.super.a1 <- try(survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]==1), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                     newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE))
    }
    estimation.Sj.a1 <- estimation.Sj.super.a1$event.SL.predict
    
    dml.train <- dml.estimation <- list("dml.rmtlj.ate.ueif"=NA, "dml.rmtlj.ate.trunc.ueif"=NA, "dml.rmtlj.ato.ueif"=NA, "dml.rmtlj.atm.ueif"=NA, "dml.rmtlj.aten.ueif"=NA)
    
    ## train
    dml.train$dml.rmtlj.ate.ueif <- double_ml_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sj.a1=train.Sj.a1)
    dml.train$dml.rmtlj.ate.trunc.ueif <- double_ml_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sj.a1=train.Sj.a1)
    dml.train$dml.rmtlj.ato.ueif <- double_ml_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sj.a1=train.Sj.a1)
    dml.train$dml.rmtlj.atm.ueif <- double_ml_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sj.a1=train.Sj.a1)
    dml.train$dml.rmtlj.aten.ueif <- double_ml_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sj.a1=train.Sj.a1)
    
    ## estimation
    dml.estimation$dml.rmtlj.ate.ueif <- double_ml_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sj.a1=estimation.Sj.a1)
    dml.estimation$dml.rmtlj.ate.trunc.ueif <- double_ml_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sj.a1=estimation.Sj.a1)
    dml.estimation$dml.rmtlj.ato.ueif <- double_ml_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sj.a1=estimation.Sj.a1)
    dml.estimation$dml.rmtlj.atm.ueif <- double_ml_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sj.a1=estimation.Sj.a1)
    dml.estimation$dml.rmtlj.aten.ueif <- double_ml_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sj.a1=estimation.Sj.a1)
    
    dml.result <- rep(list(data.frame(estimate.a0=rep(NA, 40), estimate.a1=rep(NA, 40), estimate.diff=rep(NA, 40), estimate.a0.se=rep(NA, 40), estimate.a1.se=rep(NA, 40), estimate.diff.se=rep(NA, 40))), length(taus))
    names(dml.result) <- paste0("tau", taus)
    dml.result <- lapply(dml.result, "colnames<-", paste0(colnames(dml.result[[1]]), ".v", i))
    
    dml.result <- lapply(dml.result, "rownames<-", c("dml.rmtlj.ate", "dml.rmtlj.ate.trunc", "dml.rmtlj.ato", "dml.rmtlj.atm", "dml.rmtlj.aten"))
    for (j in 1:length(taus))
    {
      #
      dml.result[[paste0("tau", taus[j])]]["dml.rmtlj.ate",] <- rmtlj_extract_avarage(train.ueif.a1=dml.train$dr.integral.mc.ht.driptw.ueif$ueif.a1, train.ueif.a0=dml.train$dr.integral.mc.ht.driptw.ueif$ueif.a0,
                                                                                      estimation.ueif.a1=dml.estimation$dr.integral.mc.ht.driptw.ueif$ueif.a1, estimation.ueif.a0=dml.estimation$dr.integral.mc.ht.driptw.ueif$ueif.a0,
                                                                                      train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
      dml.result[[paste0("tau", taus[j])]]["dml.rmtlj.ate.trunc",] <- rmtlj_extract_avarage(train.ueif.a1=dml.train$dr.integral.mc.ht.driptw.trunc.ueif$ueif.a1, train.ueif.a0=dml.train$dr.integral.mc.ht.driptw.trunc.ueif$ueif.a0,
                                                                                            estimation.ueif.a1=dml.estimation$dr.integral.mc.ht.driptw.trunc.ueif$ueif.a1, estimation.ueif.a0=dml.estimation$dr.integral.mc.ht.driptw.trunc.ueif$ueif.a0,
                                                                                            train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
      dml.result[[paste0("tau", taus[j])]]["dml.rmtlj.ato",] <- rmtlj_extract_avarage(train.ueif.a1=dml.train$dr.integral.mc.ht.drow.ueif$ueif.a1, train.ueif.a0=dml.train$dr.integral.mc.ht.drow.ueif$ueif.a0,
                                                                                      estimation.ueif.a1=dml.estimation$dr.integral.mc.ht.drow.ueif$ueif.a1, estimation.ueif.a0=dml.estimation$dr.integral.mc.ht.drow.ueif$ueif.a0,
                                                                                      train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
      dml.result[[paste0("tau", taus[j])]]["dml.rmtlj.atm",] <- rmtlj_extract_avarage(train.ueif.a1=dml.train$dr.integral.mc.ht.drmw.ueif$ueif.a1, train.ueif.a0=dml.train$dr.integral.mc.ht.drmw.ueif$ueif.a0,
                                                                                      estimation.ueif.a1=dml.estimation$dr.integral.mc.ht.drmw.ueif$ueif.a1, estimation.ueif.a0=dml.estimation$dr.integral.mc.ht.drmw.ueif$ueif.a0,
                                                                                      train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
      dml.result[[paste0("tau", taus[j])]]["dml.rmtlj.aten",] <- rmtlj_extract_avarage(train.ueif.a1=dml.train$dr.integral.mc.ht.drenw.ueif$ueif.a1, train.ueif.a0=dml.train$dr.integral.mc.ht.drenw.ueif$ueif.a0,
                                                                                       estimation.ueif.a1=dml.estimation$dr.integral.mc.ht.drenw.ueif$ueif.a1, estimation.ueif.a0=dml.estimation$dr.integral.mc.ht.drenw.ueif$ueif.a0,
                                                                                       train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    }
    # if (i%%10==0){
    print(i)
    # flush.console()
    # }
    dml.result
  }
}


# baseline
extract_linear_interp_point_estimate <- function(estimate.df, tau){
  estimate.a0 <- approx(estimate.df$s, estimate.df$estimate.a0, xout=tau)$y
  estimate.a1 <- approx(estimate.df$s, estimate.df$estimate.a1, xout=tau)$y
  estimate.diff <- estimate.a1-estimate.a0
  return(data.frame(estimate.a0, estimate.a1, estimate.diff))
}

run_rmtlj_baseline_simulation <- function(nsim, n, admin.cens, time.grid, taus, nboot) {
  baseline.result <- data.frame(matrix(NA,nrow=14,ncol=0))
  
  start.time <- Sys.time()
  for(i in 1:nsim) {
    sim.df <- dgp_competing(2*n, admin.cens, time.grid)
    
    while(max(sim.df$end.time[sim.df$a==0 & sim.df$event==0])<taus[4] | max(sim.df$end.time[sim.df$a==1 & sim.df$event==0])<taus[4] |
          max(sim.df$end.time[sim.df$a==0 & sim.df$event==1])<taus[4] | max(sim.df$end.time[sim.df$a==1 & sim.df$event==1])<taus[4]) {
      sim.df <- dgp_competing(2*n, admin.cens, time.grid)}
    
    sim.df$naive <- 1
    sim.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict(glm(a~x1+x2+x3+x4+x5+x6, data=sim.df, family=binomial()), type="response", newdata=sim.df)) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
    sim.df$iptw <- with(sim.df,1/(ps*a+(1-ps)*(1-a)))
    sim.df$iptw.trunc <- with(sim.df,ifelse(iptw>quantile(iptw,0.99),quantile(iptw,0.99),iptw))
    sim.df$ow <- with(sim.df,ifelse(a,1-ps,ps))
    sim.df$mw <- with(sim.df,pmin(ps,1-ps)*ifelse(a,1/ps,1/(1-ps)))
    sim.df$enw <- with(sim.df,-(ps*log(ps)+(1-ps)*log(1-ps))*ifelse(a,1/ps,1/(1-ps)))
    sim.df$ow.tilt <- with(sim.df, ps*(1-ps))
    sim.df$mw.tilt <- with(sim.df, pmin(ps, 1-ps))
    sim.df$enw.tilt <- with(sim.df,-(ps*log(ps)+(1-ps)*log(1-ps)))
    
    # conditional survival
    S.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event!=0))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==0,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
    S.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event!=0))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==1,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
    
    # conditional censoring
    G.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==0))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==0,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
    G.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==0))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==1,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
    
    # conditional competing
    Sj.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==0,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
    Sj.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==1,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
    
    baseline <- list("or.rmtlj.ate"=NA, "or.rmtlj.ato"=NA, "or.rmtlj.atm"=NA, "or.rmtlj.aten"=NA,
                  "ipcw.rmtlj.ate"=NA, "ipcw.rmtlj.ate.trunc"=NA, "ipcw.rmtlj.ato"=NA, "ipcw.rmtlj.atm"=NA, "ipcw.rmtlj.aten"=NA,
                  "dr.rmtlj.ate"=NA, "dr.rmtlj.ate.trunc"=NA, "dr.rmtlj.ato"=NA, "dr.rmtlj.atm"=NA, "dr.rmtlj.aten"=NA)
    
    baseline$or.rmtlj.ate <- or_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, tilt=sim.df$naive, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    baseline$or.rmtlj.ato <- or_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, tilt=sim.df$ow.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    baseline$or.rmtlj.atm <- or_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, tilt=sim.df$mw.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    baseline$or.rmtlj.aten <- or_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, tilt=sim.df$enw.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    
    baseline$ipcw.rmtlj.ate <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$iptw, G.a0=G.a0, G.a1=G.a1)
    baseline$ipcw.rmtlj.ate.trunc <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$iptw.trunc, G.a0=G.a0, G.a1=G.a1)
    baseline$ipcw.rmtlj.ato <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$ow, G.a0=G.a0, G.a1=G.a1)
    baseline$ipcw.rmtlj.atm <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$mw, G.a0=G.a0, G.a1=G.a1)
    baseline$ipcw.rmtlj.aten <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$enw, G.a0=G.a0, G.a1=G.a1)
    
    baseline$dr.rmtlj.ate <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$iptw, tilt=sim.df$naive, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    baseline$dr.rmtlj.ate.trunc <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$iptw.trunc, tilt=sim.df$naive, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    baseline$dr.rmtlj.ato <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$ow, tilt=sim.df$ow.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    baseline$dr.rmtlj.atm <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$mw, tilt=sim.df$mw.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    baseline$dr.rmtlj.aten <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$enw, tilt=sim.df$enw.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    
    baseline.point <- data.frame(matrix(NA, nrow=length(baseline), ncol=3))
    colnames(baseline.point) <- c("estimate.a0", "estimate.a1", "estimate.diff")
    rownames(baseline.point) <- names(baseline)
    
    baseline.point["or.rmtlj.ate", ] <- extract_linear_interp_point_estimate(baseline$or.rmtlj.ate, taus[4])
    baseline.point["or.rmtlj.ato", ] <- extract_linear_interp_point_estimate(baseline$or.rmtlj.ato, taus[4])
    baseline.point["or.rmtlj.atm", ] <- extract_linear_interp_point_estimate(baseline$or.rmtlj.atm, taus[4])
    baseline.point["or.rmtlj.aten", ] <- extract_linear_interp_point_estimate(baseline$or.rmtlj.aten, taus[4])
    
    baseline.point["ipcw.rmtlj.ate", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.ate, taus[4])
    baseline.point["ipcw.rmtlj.ate.trunc", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.ate.trunc, taus[4])
    baseline.point["ipcw.rmtlj.ato", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.ato, taus[4])
    baseline.point["ipcw.rmtlj.atm", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.atm, taus[4])
    baseline.point["ipcw.rmtlj.aten", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.aten, taus[4])
    
    baseline.point["dr.rmtlj.ate", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.ate, taus[4])
    baseline.point["dr.rmtlj.ate.trunc", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.ate.trunc, taus[4])
    baseline.point["dr.rmtlj.ato", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.ato, taus[4])
    baseline.point["dr.rmtlj.atm", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.atm, taus[4])
    baseline.point["dr.rmtlj.aten", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.aten, taus[4])
    
    baseline.boot <- foreach(j=1:nboot,.combine=cbind) %dopar% {
      sample.df <- sim.df[sample(1:nrow(sim.df), size = nrow(sim.df), replace = TRUE),]
      
      while(max(sample.df$end.time[sample.df$a==0 & sample.df$event==0])<taus[4] | max(sample.df$end.time[sample.df$a==1 & sample.df$event==0])<taus[4] |
            max(sample.df$end.time[sample.df$a==0 & sample.df$event==1])<taus[4] | max(sample.df$end.time[sample.df$a==1 & sample.df$event==1])<taus[4]) {
        sample.df <- sim.df[sample(1:nrow(sim.df), size = nrow(sim.df), replace = TRUE),]}
      
      sample.df$naive <- 1
      sample.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict(glm(a~x1+x2+x3+x4+x5+x6, data=sample.df, family=binomial()), type="response", newdata=sample.df)) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
      sample.df$iptw <- with(sample.df,1/(ps*a+(1-ps)*(1-a)))
      sample.df$iptw.trunc <- with(sample.df,ifelse(iptw>quantile(iptw,0.99),quantile(iptw,0.99),iptw))
      sample.df$ow <- with(sample.df,ifelse(a,1-ps,ps))
      sample.df$mw <- with(sample.df,pmin(ps,1-ps)*ifelse(a,1/ps,1/(1-ps)))
      sample.df$enw <- with(sample.df,-(ps*log(ps)+(1-ps)*log(1-ps))*ifelse(a,1/ps,1/(1-ps)))
      sample.df$ow.tilt <- with(sample.df, ps*(1-ps))
      sample.df$mw.tilt <- with(sample.df, pmin(ps, 1-ps))
      sample.df$enw.tilt <- with(sample.df,-(ps*log(ps)+(1-ps)*log(1-ps)))
      
      # conditional survival
      S.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==0,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
      S.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==1,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
      
      # conditional censoring
      G.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==0))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==0,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
      G.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==0))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==1,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
      
      # conditional competing
      Sj.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==0,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
      Sj.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==1,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
      
      baseline <- list("or.rmtlj.ate"=NA, "or.rmtlj.ato"=NA, "or.rmtlj.atm"=NA, "or.rmtlj.aten"=NA,
                       "ipcw.rmtlj.ate"=NA, "ipcw.rmtlj.ate.trunc"=NA, "ipcw.rmtlj.ato"=NA, "ipcw.rmtlj.atm"=NA, "ipcw.rmtlj.aten"=NA,
                       "dr.rmtlj.ate"=NA, "dr.rmtlj.ate.trunc"=NA, "dr.rmtlj.ato"=NA, "dr.rmtlj.atm"=NA, "dr.rmtlj.aten"=NA)
      
      baseline$or.rmtlj.ate <- or_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, tilt=sample.df$naive, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
      baseline$or.rmtlj.ato <- or_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, tilt=sample.df$ow.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
      baseline$or.rmtlj.atm <- or_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, tilt=sample.df$mw.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
      baseline$or.rmtlj.aten <- or_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, tilt=sample.df$enw.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
      
      baseline$ipcw.rmtlj.ate <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$iptw, G.a0=G.a0, G.a1=G.a1)
      baseline$ipcw.rmtlj.ate.trunc <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$iptw.trunc, G.a0=G.a0, G.a1=G.a1)
      baseline$ipcw.rmtlj.ato <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$ow, G.a0=G.a0, G.a1=G.a1)
      baseline$ipcw.rmtlj.atm <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$mw, G.a0=G.a0, G.a1=G.a1)
      baseline$ipcw.rmtlj.aten <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$enw, G.a0=G.a0, G.a1=G.a1)
      
      baseline$dr.rmtlj.ate <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$iptw, tilt=sample.df$naive, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
      baseline$dr.rmtlj.ate.trunc <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$iptw.trunc, tilt=sample.df$naive, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
      baseline$dr.rmtlj.ato <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$ow, tilt=sample.df$ow.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
      baseline$dr.rmtlj.atm <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$mw, tilt=sample.df$mw.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
      baseline$dr.rmtlj.aten <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$enw, tilt=sample.df$enw.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
      
      baseline.point <- data.frame(matrix(NA, nrow=length(baseline), ncol=3))
      colnames(baseline.point) <- c(paste0("estimate.a0.v",j), paste0("estimate.a1.v",j), paste0("estimate.diff.v",j))
      rownames(baseline.point) <- names(baseline)
      
      baseline.point["or.rmtlj.ate", ] <- extract_linear_interp_point_estimate(baseline$or.rmtlj.ate, taus[4])
      baseline.point["or.rmtlj.ato", ] <- extract_linear_interp_point_estimate(baseline$or.rmtlj.ato, taus[4])
      baseline.point["or.rmtlj.atm", ] <- extract_linear_interp_point_estimate(baseline$or.rmtlj.atm, taus[4])
      baseline.point["or.rmtlj.aten", ] <- extract_linear_interp_point_estimate(baseline$or.rmtlj.aten, taus[4])
      
      baseline.point["ipcw.rmtlj.ate", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.ate, taus[4])
      baseline.point["ipcw.rmtlj.ate.trunc", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.ate.trunc, taus[4])
      baseline.point["ipcw.rmtlj.ato", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.ato, taus[4])
      baseline.point["ipcw.rmtlj.atm", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.atm, taus[4])
      baseline.point["ipcw.rmtlj.aten", ] <- extract_linear_interp_point_estimate(baseline$ipcw.rmtlj.aten, taus[4])
      
      baseline.point["dr.rmtlj.ate", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.ate, taus[4])
      baseline.point["dr.rmtlj.ate.trunc", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.ate.trunc, taus[4])
      baseline.point["dr.rmtlj.ato", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.ato, taus[4])
      baseline.point["dr.rmtlj.atm", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.atm, taus[4])
      baseline.point["dr.rmtlj.aten", ] <- extract_linear_interp_point_estimate(baseline$dr.rmtlj.aten, taus[4])
      
      baseline.point
    }
    baseline.point$estimate.a0.se <- apply(baseline.boot[, grep("estimate.a0",names(baseline.boot))],1, sd)
    baseline.point$estimate.a1.se <- apply(baseline.boot[, grep("estimate.a1",names(baseline.boot))],1, sd)
    baseline.point$estimate.diff.se <- apply(baseline.boot[, grep("estimate.diff",names(baseline.boot))],1, sd)
    
    names(baseline.point) <- paste0(names(baseline.point),".v",i)
    baseline.result <- cbind(baseline.result, baseline.point)
    
    # if (i%%10==0){
    print(i)
    # flush.console()
    # }
  }
}


# true value
rmtlj_find_true_value <- function(true.a0, true.a1, tau)
{
  true.value <- data.frame(tau, rmtlj.a0=NA, rmtlj.ow.a0=NA, rmtlj.mw.a0=NA, rmtlj.enw.a0=NA,
                           rmtlj.a1=NA, rmtlj.ow.a1=NA, rmtlj.mw.a1=NA, rmtlj.enw.a1=NA,
                           rmtlj.diff=NA, rmtlj.ow.diff=NA, rmtlj.mw.diff=NA, rmtlj.enw.diff=NA)
  # a=0
  ind.a0 <- max(which(true.a0$s<=tau))
  true.value$rmtlj.a0 <- as.numeric(approx(x=true.a0$s[ind.a0:(ind.a0+1)], y=true.a0$rmtlj.a0[ind.a0:(ind.a0+1)], xout=tau))[2]
  true.value$rmtlj.ow.a0 <- as.numeric(approx(x=true.a0$s[ind.a0:(ind.a0+1)], y=true.a0$rmtlj.ow.a0[ind.a0:(ind.a0+1)], xout=tau))[2]
  true.value$rmtlj.mw.a0 <- as.numeric(approx(x=true.a0$s[ind.a0:(ind.a0+1)], y=true.a0$rmtlj.mw.a0[ind.a0:(ind.a0+1)], xout=tau))[2]
  true.value$rmtlj.enw.a0 <- as.numeric(approx(x=true.a0$s[ind.a0:(ind.a0+1)], y=true.a0$rmtlj.enw.a0[ind.a0:(ind.a0+1)], xout=tau))[2]
  
  # a=1
  ind.a1 <- max(which(true.a1$s<=tau))
  true.value$rmtlj.a1 <- as.numeric(approx(x=true.a1$s[ind.a1:(ind.a1+1)], y=true.a1$rmtlj.a1[ind.a1:(ind.a1+1)], xout=tau))[2]
  true.value$rmtlj.ow.a1 <- as.numeric(approx(x=true.a1$s[ind.a1:(ind.a1+1)], y=true.a1$rmtlj.ow.a1[ind.a1:(ind.a1+1)], xout=tau))[2]
  true.value$rmtlj.mw.a1 <- as.numeric(approx(x=true.a1$s[ind.a1:(ind.a1+1)], y=true.a1$rmtlj.mw.a1[ind.a1:(ind.a1+1)], xout=tau))[2]
  true.value$rmtlj.enw.a1 <- as.numeric(approx(x=true.a1$s[ind.a1:(ind.a1+1)], y=true.a1$rmtlj.enw.a1[ind.a1:(ind.a1+1)], xout=tau))[2]
  
  # diff
  true.value$rmtlj.diff <- true.value$rmtlj.a1-true.value$rmtlj.a0
  true.value$rmtlj.ow.diff <- true.value$rmtlj.ow.a1-true.value$rmtlj.ow.a0
  true.value$rmtlj.mw.diff <- true.value$rmtlj.mw.a1-true.value$rmtlj.mw.a0
  true.value$rmtlj.enw.diff <- true.value$rmtlj.enw.a1-true.value$rmtlj.enw.a0
  return(true.value)
}

rmtlj_compute_true_value <- function(time, event, tilt)
{
  n <- length(time)
  cause <- 1
  s <- sort(unique(time))
  ns <- length(s)
  ds <- diff(c(0, s))
  return(1/mean(tilt)*colMeans(matrix(tilt, ncol=ns, nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) (s[u]-pmin(time, s[u]))*(event == cause)))))
}

