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
data.df$event <- data.df$rel+2*data.df$trm # table(data.df$rel+2*data.df$trm)
# data.df$rel is the event of interest, data.df$trm is the competing event

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
train.Sj.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]==1), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                     newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
train.Sj.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]==1), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                     newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)

train.Sjbar.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]>1), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                        newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
train.Sjbar.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]>1), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                        newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)

train.all.survival.super.a0 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]!=0), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
train.G.a0 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a0$cens.SL.predict)
train.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a0$event.SL.predict)

train.all.survival.super.a1 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]!=0), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
train.G.a1 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a1$cens.SL.predict)
train.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, train.all.survival.super.a1$event.SL.predict)

# estimation
estimation.Sj.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]==1), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                          newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
estimation.Sj.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]==1), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                          newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)

estimation.Sjbar.a0 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]>1), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                             newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
estimation.Sjbar.a1 <- Winsorize(minval = 1e-3, maxval = 1, survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]>1), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                                             newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)$event.SL.predict)
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

### plugin
dml.train <- dml.estimation <- list("dml.cifj.sep.direct.astar1.ate.ueif"=NA, "dml.rmtlj.sep.direct.astar1.ate.ueif"=NA, "dml.cifj.sep.indirect.astar1.ate.ueif"=NA, "dml.rmtlj.sep.indirect.astar1.ate.ueif"=NA, "dml.cifj.sep.direct.astar0.ate.ueif"=NA, "dml.rmtlj.sep.direct.astar0.ate.ueif"=NA, "dml.cifj.sep.indirect.astar0.ate.ueif"=NA, "dml.rmtlj.sep.indirect.astar0.ate.ueif"=NA, 
                                    "dml.cifj.sep.direct.astar1.ato.ueif"=NA, "dml.rmtlj.sep.direct.astar1.ato.ueif"=NA, "dml.cifj.sep.indirect.astar1.ato.ueif"=NA, "dml.rmtlj.sep.indirect.astar1.ato.ueif"=NA, "dml.cifj.sep.direct.astar0.ato.ueif"=NA, "dml.rmtlj.sep.direct.astar0.ato.ueif"=NA, "dml.cifj.sep.indirect.astar0.ato.ueif"=NA, "dml.rmtlj.sep.indirect.astar0.ato.ueif"=NA)

## train
# ate
dml.train$dml.cifj.sep.direct.astar1.ate.ueif <- double_ml_cifj_sep_direct_astar1(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive,
                                                                               G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
dml.train$dml.rmtlj.sep.direct.astar1.ate.ueif <- double_ml_rmtlj_sep_direct_astar1(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive,
                                                                   G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)

dml.train$dml.cifj.sep.indirect.astar1.ate.ueif <- double_ml_cifj_sep_indirect_astar1(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive,
                                                                                   G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
dml.train$dml.rmtlj.sep.indirect.astar1.ate.ueif <- double_ml_rmtlj_sep_indirect_astar1(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive,
                                                                       G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)

dml.train$dml.cifj.sep.direct.astar0.ate.ueif <- double_ml_cifj_sep_direct_astar0(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive,
                                                                               G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
dml.train$dml.rmtlj.sep.direct.astar0.ate.ueif <- double_ml_rmtlj_sep_direct_astar0(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive,
                                                                            G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)

dml.train$dml.cifj.sep.indirect.astar0.ate.ueif <- double_ml_cifj_sep_indirect_astar0(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive,
                                                                                   G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
dml.train$dml.rmtlj.sep.indirect.astar0.ate.ueif <- double_ml_rmtlj_sep_indirect_astar0(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive,
                                                                                G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)

## train
# ato
dml.train$dml.cifj.sep.direct.astar1.ato.ueif <- double_ml_cifj_sep_direct_astar1(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt,
                                                                             G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
dml.train$dml.rmtlj.sep.direct.astar1.ato.ueif <- double_ml_rmtlj_sep_direct_astar1(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt,
                                                                          G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)

dml.train$dml.cifj.sep.indirect.astar1.ato.ueif <- double_ml_cifj_sep_indirect_astar1(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt,
                                                                                 G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
dml.train$dml.rmtlj.sep.indirect.astar1.ato.ueif <- double_ml_rmtlj_sep_indirect_astar1(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt,
                                                                              G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)

dml.train$dml.cifj.sep.direct.astar0.ato.ueif <- double_ml_cifj_sep_direct_astar0(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt,
                                                                             G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
dml.train$dml.rmtlj.sep.direct.astar0.ato.ueif <- double_ml_rmtlj_sep_direct_astar0(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt,
                                                                          G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)

dml.train$dml.cifj.sep.indirect.astar0.ato.ueif <- double_ml_cifj_sep_indirect_astar0(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt,
                                                                                 G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)
dml.train$dml.rmtlj.sep.indirect.astar0.ato.ueif <- double_ml_rmtlj_sep_indirect_astar0(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt,
                                                                              G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1, Sj.a0=train.Sj.a0, Sjbar.a0=train.Sjbar.a0, Sj.a1=train.Sj.a1, Sjbar.a1=train.Sjbar.a1, freq.time=NULL, admin.cens)

## estimation
# ate
dml.estimation$dml.cifj.sep.direct.astar1.ate.ueif <- double_ml_cifj_sep_direct_astar1(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive,
                                                                                    G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
dml.estimation$dml.rmtlj.sep.direct.astar1.ate.ueif <- double_ml_rmtlj_sep_direct_astar1(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive,
                                                                                 G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)

dml.estimation$dml.cifj.sep.indirect.astar1.ate.ueif <- double_ml_cifj_sep_indirect_astar1(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive,
                                                                                        G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
dml.estimation$dml.rmtlj.sep.indirect.astar1.ate.ueif <- double_ml_rmtlj_sep_indirect_astar1(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive,
                                                                                     G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)

dml.estimation$dml.cifj.sep.direct.astar0.ate.ueif <- double_ml_cifj_sep_direct_astar0(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive,
                                                                                    G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
dml.estimation$dml.rmtlj.sep.direct.astar0.ate.ueif <- double_ml_rmtlj_sep_direct_astar0(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive,
                                                                                 G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)

dml.estimation$dml.cifj.sep.indirect.astar0.ate.ueif <- double_ml_cifj_sep_indirect_astar0(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive,
                                                                                        G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
dml.estimation$dml.rmtlj.sep.indirect.astar0.ate.ueif <- double_ml_rmtlj_sep_indirect_astar0(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive,
                                                                                     G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)

## estimation
# ato
dml.estimation$dml.cifj.sep.direct.astar1.ato.ueif <- double_ml_cifj_sep_direct_astar1(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt,
                                                                                  G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
dml.estimation$dml.rmtlj.sep.direct.astar1.ato.ueif <- double_ml_rmtlj_sep_direct_astar1(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt,
                                                                               G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)

dml.estimation$dml.cifj.sep.indirect.astar1.ato.ueif <- double_ml_cifj_sep_indirect_astar1(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt,
                                                                                      G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
dml.estimation$dml.rmtlj.sep.indirect.astar1.ato.ueif <- double_ml_rmtlj_sep_indirect_astar1(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt,
                                                                                   G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)

dml.estimation$dml.cifj.sep.direct.astar0.ato.ueif <- double_ml_cifj_sep_direct_astar0(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt,
                                                                                  G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
dml.estimation$dml.rmtlj.sep.direct.astar0.ato.ueif <- double_ml_rmtlj_sep_direct_astar0(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt,
                                                                               G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)

dml.estimation$dml.cifj.sep.indirect.astar0.ato.ueif <- double_ml_cifj_sep_indirect_astar0(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt,
                                                                                      G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)
dml.estimation$dml.rmtlj.sep.indirect.astar0.ato.ueif <- double_ml_rmtlj_sep_indirect_astar0(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt,
                                                                                   G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1, Sj.a0=estimation.Sj.a0, Sjbar.a0=estimation.Sjbar.a0, Sj.a1=estimation.Sj.a1, Sjbar.a1=estimation.Sjbar.a1, freq.time=NULL, admin.cens)

rm(train.G.a0, train.G.a1, train.S.a0, train.S.a1, train.Sj.a0, train.Sj.a1, train.Sjbar.a0, train.Sjbar.a1, train.all.survival.super.a0, train.all.survival.super.a1,
   estimation.G.a0, estimation.G.a1, estimation.S.a0, estimation.S.a1, estimation.Sj.a0, estimation.Sj.a1, estimation.Sjbar.a0, estimation.Sjbar.a1, estimation.all.survival.super.a0, estimation.all.survival.super.a1)
gc()

# inference
start.time <- Sys.time()

relapse.separable.ls <- list("dml.cifj.sep.direct.astar1.ate"=NA, "dml.rmtlj.sep.direct.astar1.ate"=NA, "dml.cifj.sep.indirect.astar1.ate"=NA, "dml.rmtlj.sep.indirect.astar1.ate"=NA, "dml.cifj.sep.direct.astar0.ate"=NA, "dml.rmtlj.sep.direct.astar0.ate"=NA, "dml.cifj.sep.indirect.astar0.ate"=NA, "dml.rmtlj.sep.indirect.astar0.ate"=NA,
                             "dml.cifj.sep.direct.astar1.ato"=NA, "dml.rmtlj.sep.direct.astar1.ato"=NA, "dml.cifj.sep.indirect.astar1.ato"=NA, "dml.rmtlj.sep.indirect.astar1.ato"=NA, "dml.cifj.sep.direct.astar0.ato"=NA, "dml.rmtlj.sep.direct.astar0.ato"=NA, "dml.cifj.sep.indirect.astar0.ato"=NA, "dml.rmtlj.sep.indirect.astar0.ato"=NA)

# cifj, ate
relapse.separable.ls$dml.cifj.sep.direct.astar1.ate <- cifj_inference(npath=npath, train.ueif.a1=f.train$sep.direct.astar1.cifj.driptw.ueif$ueif.a1, train.ueif.a0=f.train$sep.direct.astar1.cifj.driptw.ueif$ueif.a0,
                                                                                         estimation.ueif.a1=f.estimation$sep.direct.astar1.cifj.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.direct.astar1.cifj.driptw.ueif$ueif.a0,
                                                                                         train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.cifj.sep.indirect.astar1.ate <- cifj_inference(npath=npath, train.ueif.a1=f.train$sep.indirect.astar1.cifj.driptw.ueif$ueif.a1, train.ueif.a0=f.train$sep.indirect.astar1.cifj.driptw.ueif$ueif.a0,
                                                                   estimation.ueif.a1=f.estimation$sep.indirect.astar1.cifj.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.indirect.astar1.cifj.driptw.ueif$ueif.a0,
                                                                   train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.cifj.sep.direct.astar0.ate <- cifj_inference(npath=npath, train.ueif.a1=f.train$sep.direct.astar0.cifj.driptw.ueif$ueif.a1, train.ueif.a0=f.train$sep.direct.astar0.cifj.driptw.ueif$ueif.a0,
                                                                 estimation.ueif.a1=f.estimation$sep.direct.astar0.cifj.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.direct.astar0.cifj.driptw.ueif$ueif.a0,
                                                                 train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.cifj.sep.indirect.astar0.ate <- cifj_inference(npath=npath, train.ueif.a1=f.train$sep.indirect.astar0.cifj.driptw.ueif$ueif.a1, train.ueif.a0=f.train$sep.indirect.astar0.cifj.driptw.ueif$ueif.a0,
                                                                   estimation.ueif.a1=f.estimation$sep.indirect.astar0.cifj.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.indirect.astar0.cifj.driptw.ueif$ueif.a0,
                                                                   train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))

# rmtlj, ate
relapse.separable.ls$dml.rmtlj.sep.direct.astar1.ate <- rmtlj_inference(npath=npath, train.ueif.a1=f.train$sep.direct.astar1.rmtlj.driptw.ueif$ueif.a1, train.ueif.a0=f.train$sep.direct.astar1.rmtlj.driptw.ueif$ueif.a0,
                                                                estimation.ueif.a1=f.estimation$sep.direct.astar1.rmtlj.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.direct.astar1.rmtlj.driptw.ueif$ueif.a0,
                                                                train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.rmtlj.sep.indirect.astar1.ate <- rmtlj_inference(npath=npath, train.ueif.a1=f.train$sep.indirect.astar1.rmtlj.driptw.ueif$ueif.a1, train.ueif.a0=f.train$sep.indirect.astar1.rmtlj.driptw.ueif$ueif.a0,
                                                                  estimation.ueif.a1=f.estimation$sep.indirect.astar1.rmtlj.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.indirect.astar1.rmtlj.driptw.ueif$ueif.a0,
                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.rmtlj.sep.direct.astar0.ate <- rmtlj_inference(npath=npath, train.ueif.a1=f.train$sep.direct.astar0.rmtlj.driptw.ueif$ueif.a1, train.ueif.a0=f.train$sep.direct.astar0.rmtlj.driptw.ueif$ueif.a0,
                                                                   estimation.ueif.a1=f.estimation$sep.direct.astar0.rmtlj.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.direct.astar0.rmtlj.driptw.ueif$ueif.a0,
                                                                   train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.rmtlj.sep.indirect.astar0.ate <- rmtlj_inference(npath=npath, train.ueif.a1=f.train$sep.indirect.astar0.rmtlj.driptw.ueif$ueif.a1, train.ueif.a0=f.train$sep.indirect.astar0.rmtlj.driptw.ueif$ueif.a0,
                                                                     estimation.ueif.a1=f.estimation$sep.indirect.astar0.rmtlj.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.indirect.astar0.rmtlj.driptw.ueif$ueif.a0,
                                                                     train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))

# cifj, ato
relapse.separable.ls$dml.cifj.sep.direct.astar1.ato <- cifj_inference(npath=npath, train.ueif.a1=f.train$sep.direct.astar1.cifj.drow.ueif$ueif.a1, train.ueif.a0=f.train$sep.direct.astar1.cifj.drow.ueif$ueif.a0,
                                                               estimation.ueif.a1=f.estimation$sep.direct.astar1.cifj.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.direct.astar1.cifj.drow.ueif$ueif.a0,
                                                               train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.cifj.sep.indirect.astar1.ato <- cifj_inference(npath=npath, train.ueif.a1=f.train$sep.indirect.astar1.cifj.drow.ueif$ueif.a1, train.ueif.a0=f.train$sep.indirect.astar1.cifj.drow.ueif$ueif.a0,
                                                                 estimation.ueif.a1=f.estimation$sep.indirect.astar1.cifj.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.indirect.astar1.cifj.drow.ueif$ueif.a0,
                                                                 train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.cifj.sep.direct.astar0.ato <- cifj_inference(npath=npath, train.ueif.a1=f.train$sep.direct.astar0.cifj.drow.ueif$ueif.a1, train.ueif.a0=f.train$sep.direct.astar0.cifj.drow.ueif$ueif.a0,
                                                               estimation.ueif.a1=f.estimation$sep.direct.astar0.cifj.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.direct.astar0.cifj.drow.ueif$ueif.a0,
                                                               train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.cifj.sep.indirect.astar0.ato <- cifj_inference(npath=npath, train.ueif.a1=f.train$sep.indirect.astar0.cifj.drow.ueif$ueif.a1, train.ueif.a0=f.train$sep.indirect.astar0.cifj.drow.ueif$ueif.a0,
                                                                 estimation.ueif.a1=f.estimation$sep.indirect.astar0.cifj.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.indirect.astar0.cifj.drow.ueif$ueif.a0,
                                                                 train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))

# rmtlj, ato
relapse.separable.ls$dml.rmtlj.sep.direct.astar1.ato <- rmtlj_inference(npath=npath, train.ueif.a1=f.train$sep.direct.astar1.rmtlj.drow.ueif$ueif.a1, train.ueif.a0=f.train$sep.direct.astar1.rmtlj.drow.ueif$ueif.a0,
                                                                 estimation.ueif.a1=f.estimation$sep.direct.astar1.rmtlj.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.direct.astar1.rmtlj.drow.ueif$ueif.a0,
                                                                 train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.rmtlj.sep.indirect.astar1.ato <- rmtlj_inference(npath=npath, train.ueif.a1=f.train$sep.indirect.astar1.rmtlj.drow.ueif$ueif.a1, train.ueif.a0=f.train$sep.indirect.astar1.rmtlj.drow.ueif$ueif.a0,
                                                                   estimation.ueif.a1=f.estimation$sep.indirect.astar1.rmtlj.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.indirect.astar1.rmtlj.drow.ueif$ueif.a0,
                                                                   train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.rmtlj.sep.direct.astar0.ato <- rmtlj_inference(npath=npath, train.ueif.a1=f.train$sep.direct.astar0.rmtlj.drow.ueif$ueif.a1, train.ueif.a0=f.train$sep.direct.astar0.rmtlj.drow.ueif$ueif.a0,
                                                                 estimation.ueif.a1=f.estimation$sep.direct.astar0.rmtlj.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.direct.astar0.rmtlj.drow.ueif$ueif.a0,
                                                                 train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
relapse.separable.ls$dml.rmtlj.sep.indirect.astar0.ato <- rmtlj_inference(npath=npath, train.ueif.a1=f.train$sep.indirect.astar0.rmtlj.drow.ueif$ueif.a1, train.ueif.a0=f.train$sep.indirect.astar0.rmtlj.drow.ueif$ueif.a0,
                                                                   estimation.ueif.a1=f.estimation$sep.indirect.astar0.rmtlj.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$sep.indirect.astar0.rmtlj.drow.ueif$ueif.a0,
                                                                   train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))

end.time <- Sys.time()
end.time - start.time

save(relapse.separable.ls, file="relapse_separable_wate.RData")
