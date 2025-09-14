## RAM: ; core: ;  hrs; testing
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

library(grf)
library(SuperLearner)
library(survSuperLearner)
library(xgboost)
# library(bartMachine)
# library(KernelKnn)
library(earth)
library(fdrtool)
library(simecol)
library(MASS)
library(data.table)
library(ggplot2)
# library(mvnfast)
# library(lqmm)
library(Rfast)
library(haven)

registerDoParallel(cores=8)

tau <- 4*365.25
a.sl.lib <- c("SL.mean","SL.glm","SL.xgboost","SL.ranger","SL.nnet","SL.rpartPrune","SL.step","SL.glmnet","SL.earth") # ,"SL.kernelKnn","SL.gam"
event.sl.lib <- cens.sl.lib <- c("survSL.km","survSL.coxph","survSL.rfsrc","survSL.loglogreg","survSL.expreg","survSL.weibreg") #  event.sl.lib <- cens.sl.lib <- c("survSL.km","survSL.coxph") "survSL.pchreg": too slow; # ; "survSL.gam","survSL.pchSL"
weights.sl.lib <- c("SL.glm","SL.xgboost","SL.ranger","SL.nnet","SL.earth","SL.glmnet") # c("SL.bayesglm","SL.biglasso","SL.cforest","SL.earth","SL.gam","SL.glm","SL.glmnet","SL.loess","SL.nnet","SL.nnls","SL.ranger","SL.xgboost") # "SL.bartMachine","SL.gbm","SL.caret","SL.gam"
npath <- 1000

####################################################################################
data.df <- read_sas("P-5370 (LK16-02)/lk1602public.sas7bdat")
# data.df <- data.df[(data.df$graftype != 23), ]
# data.df$graftype[data.df$graftype==22] <- 0 # 0: peripheral blood; 1: bone marrow

data.df$end.time <- round(pmin(data.df$intxrel, data.df$intxsurv)*30.4167, digits = 0)
data.df <- data.df[!is.na(data.df$end.time), ]
data.df$event <- data.df$trm+data.df$rel # composite outcome of relaspe and mortality, table(data.df$rel+data.df$trm)

data.df$risk_level <- ifelse(data.df$strataf %in% c(1, 3), 0, 1) # 1 and 3 intermediate risk

# intxagvhd, intxcgvhd, gvhdgpn, agvhd: post treatment
categorical.confounders <- c("drsex","disease","distatus","amstatprgp","kps","wbcdxgp","atggrp",
                             "gvhdgpf","drcmvprf","donorgp_final","agegpff","graftype","risk_level") # discuss
# "condgrp": finer grain of treatment intensity

data.df <- data.df[data.df$drsex!=99, ]
data.df <- data.df[data.df$atggrp!=99, ]
data.df <- data.df[data.df$kps!=99, ]
# data.df <- data.df[data.df$condintb!=99, ] # condintb = 99 only exists for a=0
# data.df <- data.df[data.df$insurgrp!=99, ]
# data.df <- data.df[data.df$marstat!=3, ]
# data.df <- data.df[data.df$marstat!=5, ]
# data.df <- data.df[data.df$donorgp!=6, ]
# data.df <- data.df[data.df$cgvhd1y!=6, ]
# data.df <- data.df[data.df$cgvhd1y!=99, ]

data.df[, categorical.confounders] <- lapply(data.df[, categorical.confounders], as.factor)

confounders <- c("yeartx", categorical.confounders)

# data.df <- data.df[complete.cases(data.df[, confounders]), ]
# data.df <- data.df[1:2000,]
time.grid <- 1
admin.cens <- max(data.df$end.time)
n <- nrow(data.df)
data.df$a <- ifelse(data.df$strataf %in% c(1, 2), 0, 1) # strataf as treatment, group 1 and 2; 3 and 4
data.df$id <- 1:n

####################################################################################
## first split
train.ind <- sample(1:n, floor(n/2), replace=FALSE) # sample training set
train.df <- data.df[train.ind,]
estimation.df <- data.df[setdiff(1:n, train.ind),]

# train
start.time <- Sys.time()
train.all.survival.super.a0 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]!=0), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
train.G.a0 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a0$cens.SL.predict)
train.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a0$event.SL.predict)

train.all.survival.super.a1 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]!=0), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
train.G.a1 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a1$cens.SL.predict)
train.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a1$event.SL.predict)

# estimation
estimation.all.survival.super.a0 <- survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]!=0), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                     newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
estimation.G.a0 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a0$cens.SL.predict)
estimation.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a0$event.SL.predict)

estimation.all.survival.super.a1 <- survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]!=0), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                     newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
estimation.G.a1 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a1$cens.SL.predict)
estimation.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, estimation.all.survival.super.a1$event.SL.predict)

end.time <- Sys.time()
end.time - start.time

# plugin
dml.train <- dml.estimation <- list("dml.cifj.ate.ueif"=NA, "dml.rmtlj.ate.ueif"=NA, "dml.cifj.ato.ueif"=NA, "dml.rmtlj.ato.ueif"=NA)

dml.train$dml.cifj.ate.ueif <- double_ml_cifj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.S.a0, Sj.a1=train.S.a1)
dml.train$dml.rmtlj.ate.ueif <- double_ml_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.S.a0, Sj.a1=train.S.a1)

dml.train$dml.cifj.ato.ueif <- double_ml_cifj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.S.a0, Sj.a1=train.S.a1)
dml.train$dml.rmtlj.ato.ueif <- double_ml_rmtlj(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.S.a0, Sj.a1=train.S.a1)

dml.estimation$dml.cifj.ate.ueif <- double_ml_cifj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.S.a0, Sj.a1=estimation.S.a1)
dml.estimation$dml.rmtlj.ate.ueif <- double_ml_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.S.a0, Sj.a1=estimation.S.a1)

dml.estimation$dml.cifj.ato.ueif <- double_ml_cifj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.S.a0, Sj.a1=estimation.S.a1)
dml.estimation$dml.rmtlj.ato.ueif <- double_ml_rmtlj(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.S.a0, Sj.a1=estimation.S.a1)

rm(train.G.a0, train.G.a1, train.S.a0, train.S.a1, train.all.survival.super.a0, train.all.survival.super.a1, estimation.G.a0, estimation.G.a1, estimation.S.a0, estimation.S.a1, estimation.all.survival.super.a0, estimation.all.survival.super.a1)
gc()

# inference
start.time <- Sys.time()

composite.ls <- list("dml.cifj.ate"=NA, "dml.rmtlj.ate"=NA, "dml.cifj.ato"=NA, "dml.rmtlj.ato"=NA)

composite.ls$dml.cifj.ate <- cifj_inference(npath=npath, train.ueif.a1=dml.train$dml.cifj.ate.ueif$ueif.a1, train.ueif.a0=dml.train$dml.cifj.ate.ueif$ueif.a0,
                                            estimation.ueif.a1=dml.estimation$dml.cifj.ate.ueif$ueif.a1, estimation.ueif.a0=dml.estimation$dml.cifj.ate.ueif$ueif.a0,
                                            train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
composite.ls$dml.rmtlj.ate <- rmtlj_inference(npath=npath, train.ueif.a1=dml.train$dml.rmtlj.ate.ueif$ueif.a1, train.ueif.a0=dml.train$dml.rmtlj.ate.ueif$ueif.a0,
                                              estimation.ueif.a1=dml.estimation$dml.rmtlj.ate.ueif$ueif.a1, estimation.ueif.a0=dml.estimation$dml.rmtlj.ate.ueif$ueif.a0,
                                              train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))

composite.ls$dml.cifj.ato <- cifj_inference(npath=npath, train.ueif.a1=dml.train$dml.cifj.ato.ueif$ueif.a1, train.ueif.a0=dml.train$dml.cifj.ato.ueif$ueif.a0,
                                            estimation.ueif.a1=dml.estimation$dml.cifj.ato.ueif$ueif.a1, estimation.ueif.a0=dml.estimation$dml.cifj.ato.ueif$ueif.a0,
                                            train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
composite.ls$dml.rmtlj.ato <- rmtlj_inference(npath=npath, train.ueif.a1=dml.train$dml.rmtlj.ato.ueif$ueif.a1, train.ueif.a0=dml.train$dml.rmtlj.ato.ueif$ueif.a0,
                                              estimation.ueif.a1=dml.estimation$dml.rmtlj.ato.ueif$ueif.a1, estimation.ueif.a0=dml.estimation$dml.rmtlj.ato.ueif$ueif.a0,
                                              train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))

end.time <- Sys.time()
end.time - start.time

save(composite.ls, file="composite_wate.RData")
