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

registerDoParallel(cores=102)

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
data.df$naive <- 1

#
data.Sj.a0 <- Winsorize(minval = 1e-3, maxval = 1, matrix(summary(survfit(Surv(end.time[data.df$a==0], 1*(event[data.df$a==0]==1)) ~ 1, data=data.df), times=sort(unique(data.df$end.time)), extend=TRUE)$surv, nrow=nrow(data.df), ncol=n_distinct(data.df$end.time), byrow=TRUE))
data.Sj.a1 <- Winsorize(minval = 1e-3, maxval = 1, matrix(summary(survfit(Surv(end.time[data.df$a==1], 1*(event[data.df$a==1]==1)) ~ 1, data=data.df), times=sort(unique(data.df$end.time)), extend=TRUE)$surv, nrow=nrow(data.df), ncol=n_distinct(data.df$end.time), byrow=TRUE))

data.Sjbar.a0 <- Winsorize(minval = 1e-3, maxval = 1, matrix(summary(survfit(Surv(end.time[data.df$a==0], 1*(event[data.df$a==0]>1)) ~ 1, data=data.df), times=sort(unique(data.df$end.time)), extend=TRUE)$surv, nrow=nrow(data.df), ncol=n_distinct(data.df$end.time), byrow=TRUE))
data.Sjbar.a1 <- Winsorize(minval = 1e-3, maxval = 1, matrix(summary(survfit(Surv(end.time[data.df$a==1], 1*(event[data.df$a==1]>1)) ~ 1, data=data.df), times=sort(unique(data.df$end.time)), extend=TRUE)$surv, nrow=nrow(data.df), ncol=n_distinct(data.df$end.time), byrow=TRUE))

data.S.a0 <- Winsorize(minval = 1e-3, maxval = 1, matrix(summary(survfit(Surv(end.time[data.df$a==0], 1*(event[data.df$a==0]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time)), extend=TRUE)$surv, nrow=nrow(data.df), ncol=n_distinct(data.df$end.time), byrow=TRUE))
data.G.a0 <- Winsorize(minval = 1e-3, maxval = 1, matrix(summary(survfit(Surv(end.time[data.df$a==0], 1*(event[data.df$a==0]==0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time)), extend=TRUE)$surv, nrow=nrow(data.df), ncol=n_distinct(data.df$end.time), byrow=TRUE))

data.S.a1 <- Winsorize(minval = 1e-3, maxval = 1, matrix(summary(survfit(Surv(end.time[data.df$a==1], 1*(event[data.df$a==1]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time)), extend=TRUE)$surv, nrow=nrow(data.df), ncol=n_distinct(data.df$end.time), byrow=TRUE))
data.G.a1 <- Winsorize(minval = 1e-3, maxval = 1, matrix(summary(survfit(Surv(end.time[data.df$a==1], 1*(event[data.df$a==1]==0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time)), extend=TRUE)$surv, nrow=nrow(data.df), ncol=n_distinct(data.df$end.time), byrow=TRUE))

# plugin
crude <- crude.train <- crude.estimation <- list("crude.cifj.sep.direct.astar1.ate.ueif"=NA, "crude.rmtlj.sep.direct.astar1.ate.ueif"=NA, "crude.cifj.sep.indirect.astar1.ate.ueif"=NA, "crude.rmtlj.sep.indirect.astar1.ate.ueif"=NA,
                                                 "crude.cifj.sep.direct.astar0.ate.ueif"=NA, "crude.rmtlj.sep.direct.astar0.ate.ueif"=NA, "crude.cifj.sep.indirect.astar0.ate.ueif"=NA, "crude.rmtlj.sep.indirect.astar0.ate.ueif"=NA)

crude$crude.cifj.sep.direct.astar1.ate.ueif <- double_ml_cifj_sep_direct_astar1(id=data.df$id, a=data.df$a, time=data.df$end.time, event=data.df$event, bw=data.df$iptw, tilt=data.df$naive,
                                                                                  G.a0=data.G.a0, G.a1=data.G.a1, S.a0=data.S.a0, S.a1=data.S.a1, Sj.a0=data.Sj.a0, Sjbar.a0=data.Sjbar.a0, Sj.a1=data.Sj.a1, Sjbar.a1=data.Sjbar.a1, freq.time=NULL, admin.cens)
crude$crude.rmtlj.sep.direct.astar1.ate.ueif <- double_ml_rmtlj_sep_direct_astar1(id=data.df$id, a=data.df$a, time=data.df$end.time, event=data.df$event, bw=data.df$iptw, tilt=data.df$naive,
                                                                                    G.a0=data.G.a0, G.a1=data.G.a1, S.a0=data.S.a0, S.a1=data.S.a1, Sj.a0=data.Sj.a0, Sjbar.a0=data.Sjbar.a0, Sj.a1=data.Sj.a1, Sjbar.a1=data.Sjbar.a1, freq.time=NULL, admin.cens)

crude$crude.cifj.sep.indirect.astar1.ate.ueif <- double_ml_cifj_sep_indirect_astar1(id=data.df$id, a=data.df$a, time=data.df$end.time, event=data.df$event, bw=data.df$iptw, tilt=data.df$naive,
                                                                                      G.a0=data.G.a0, G.a1=data.G.a1, S.a0=data.S.a0, S.a1=data.S.a1, Sj.a0=data.Sj.a0, Sjbar.a0=data.Sjbar.a0, Sj.a1=data.Sj.a1, Sjbar.a1=data.Sjbar.a1, freq.time=NULL, admin.cens)
crude$crude.rmtlj.sep.indirect.astar1.ate.ueif <- double_ml_rmtlj_sep_indirect_astar1(id=data.df$id, a=data.df$a, time=data.df$end.time, event=data.df$event, bw=data.df$iptw, tilt=data.df$naive,
                                                                                        G.a0=data.G.a0, G.a1=data.G.a1, S.a0=data.S.a0, S.a1=data.S.a1, Sj.a0=data.Sj.a0, Sjbar.a0=data.Sjbar.a0, Sj.a1=data.Sj.a1, Sjbar.a1=data.Sjbar.a1, freq.time=NULL, admin.cens)

crude$crude.cifj.sep.direct.astar0.ate.ueif <- double_ml_cifj_sep_direct_astar0(id=data.df$id, a=data.df$a, time=data.df$end.time, event=data.df$event, bw=data.df$iptw, tilt=data.df$naive,
                                                                                  G.a0=data.G.a0, G.a1=data.G.a1, S.a0=data.S.a0, S.a1=data.S.a1, Sj.a0=data.Sj.a0, Sjbar.a0=data.Sjbar.a0, Sj.a1=data.Sj.a1, Sjbar.a1=data.Sjbar.a1, freq.time=NULL, admin.cens)
crude$crude.rmtlj.sep.direct.astar0.ate.ueif <- double_ml_rmtlj_sep_direct_astar0(id=data.df$id, a=data.df$a, time=data.df$end.time, event=data.df$event, bw=data.df$iptw, tilt=data.df$naive,
                                                                                    G.a0=data.G.a0, G.a1=data.G.a1, S.a0=data.S.a0, S.a1=data.S.a1, Sj.a0=data.Sj.a0, Sjbar.a0=data.Sjbar.a0, Sj.a1=data.Sj.a1, Sjbar.a1=data.Sjbar.a1, freq.time=NULL, admin.cens)

crude$crude.cifj.sep.indirect.astar0.ate.ueif <- double_ml_cifj_sep_indirect_astar0(id=data.df$id, a=data.df$a, time=data.df$end.time, event=data.df$event, bw=data.df$iptw, tilt=data.df$naive,
                                                                                      G.a0=data.G.a0, G.a1=data.G.a1, S.a0=data.S.a0, S.a1=data.S.a1, Sj.a0=data.Sj.a0, Sjbar.a0=data.Sjbar.a0, Sj.a1=data.Sj.a1, Sjbar.a1=data.Sjbar.a1, freq.time=NULL, admin.cens)
crude$crude.rmtlj.sep.indirect.astar0.ate.ueif <- double_ml_rmtlj_sep_indirect_astar0(id=data.df$id, a=data.df$a, time=data.df$end.time, event=data.df$event, bw=data.df$iptw, tilt=data.df$naive,
                                                                                        G.a0=data.G.a0, G.a1=data.G.a1, S.a0=data.S.a0, S.a1=data.S.a1, Sj.a0=data.Sj.a0, Sjbar.a0=data.Sjbar.a0, Sj.a1=data.Sj.a1, Sjbar.a1=data.Sjbar.a1, freq.time=NULL, admin.cens)

rm(data.S.a0, data.G.a0, data.Sj.a0, data.Sjbar.a0, data.S.a1, data.G.a1, data.Sj.a1, data.Sjbar.a1)
gc()

# inference
#
crude.train$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a0 <- crude$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a0[1:ceiling(n/2), ]
crude.train$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a1 <- crude$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a1[1:ceiling(n/2), ]

crude.estimation$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a0 <- crude$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a0[setdiff(1:n, 1:ceiling(n/2)), ]
crude.estimation$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a1 <- crude$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a1[setdiff(1:n, 1:ceiling(n/2)), ]

crude.train$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a0 <- crude$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a0[1:ceiling(n/2), ]
crude.train$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a1 <- crude$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a1[1:ceiling(n/2), ]

crude.estimation$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a0 <- crude$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a0[setdiff(1:n, 1:ceiling(n/2)), ]
crude.estimation$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a1 <- crude$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a1[setdiff(1:n, 1:ceiling(n/2)), ]

#
crude.train$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a0 <- crude$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a0[1:ceiling(n/2), ]
crude.train$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a1 <- crude$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a1[1:ceiling(n/2), ]

crude.estimation$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a0 <- crude$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a0[setdiff(1:n, 1:ceiling(n/2)), ]
crude.estimation$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a1 <- crude$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a1[setdiff(1:n, 1:ceiling(n/2)), ]

crude.train$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a0 <- crude$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a0[1:ceiling(n/2), ]
crude.train$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a1 <- crude$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a1[1:ceiling(n/2), ]

crude.estimation$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a0 <- crude$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a0[setdiff(1:n, 1:ceiling(n/2)), ]
crude.estimation$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a1 <- crude$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a1[setdiff(1:n, 1:ceiling(n/2)), ]

#
crude.train$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a0 <- crude$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a0[1:ceiling(n/2), ]
crude.train$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a1 <- crude$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a1[1:ceiling(n/2), ]

crude.estimation$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a0 <- crude$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a0[setdiff(1:n, 1:ceiling(n/2)), ]
crude.estimation$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a1 <- crude$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a1[setdiff(1:n, 1:ceiling(n/2)), ]

crude.train$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a0 <- crude$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a0[1:ceiling(n/2), ]
crude.train$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a1 <- crude$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a1[1:ceiling(n/2), ]

crude.estimation$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a0 <- crude$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a0[setdiff(1:n, 1:ceiling(n/2)), ]
crude.estimation$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a1 <- crude$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a1[setdiff(1:n, 1:ceiling(n/2)), ]

#
crude.train$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a0 <- crude$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a0[1:ceiling(n/2), ]
crude.train$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a1 <- crude$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a1[1:ceiling(n/2), ]

crude.estimation$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a0 <- crude$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a0[setdiff(1:n, 1:ceiling(n/2)), ]
crude.estimation$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a1 <- crude$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a1[setdiff(1:n, 1:ceiling(n/2)), ]

crude.train$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a0 <- crude$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a0[1:ceiling(n/2), ]
crude.train$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a1 <- crude$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a1[1:ceiling(n/2), ]

crude.estimation$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a0 <- crude$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a0[setdiff(1:n, 1:ceiling(n/2)), ]
crude.estimation$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a1 <- crude$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a1[setdiff(1:n, 1:ceiling(n/2)), ]

#
relapse.separable.crude.ls <- list("crude.cifj.sep.direct.astar1.ate"=NA, "crude.rmtlj.sep.direct.astar1.ate"=NA, "crude.cifj.sep.indirect.astar1.ate"=NA, "crude.rmtlj.sep.indirect.astar1.ate"=NA,
                                   "crude.cifj.sep.direct.astar0.ate"=NA, "crude.rmtlj.sep.direct.astar0.ate"=NA, "crude.cifj.sep.indirect.astar0.ate"=NA, "crude.rmtlj.sep.indirect.astar0.ate"=NA)

#
relapse.separable.crude.ls$crude.cifj.sep.direct.astar1.ate <- cifj_inference(npath=npath, train.ueif.a1=crude.train$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a1, train.ueif.a0=crude.train$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a0,
                                                                              estimation.ueif.a1=crude.estimation$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a1, estimation.ueif.a0=crude.estimation$crude.cifj.sep.direct.astar1.ate.ueif$ueif.a0,
                                                                              data.time=sort(unique(data.df$end.time)), estimation.time=sort(unique(data.df$end.time)))

relapse.separable.crude.ls$crude.rmtlj.sep.direct.astar1.ate <- rmtlj_inference(npath=npath, train.ueif.a1=crude.train$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a1, train.ueif.a0=crude.train$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a0,
                                                                                estimation.ueif.a1=crude.estimation$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a1, estimation.ueif.a0=crude.estimation$crude.rmtlj.sep.direct.astar1.ate.ueif$ueif.a0,
                                                                                data.time=sort(unique(data.df$end.time)), estimation.time=sort(unique(data.df$end.time)))
#
relapse.separable.crude.ls$crude.cifj.sep.indirect.astar1.ate <- cifj_inference(npath=npath, train.ueif.a1=crude.train$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a1, train.ueif.a0=crude.train$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a0,
                                                                                estimation.ueif.a1=crude.estimation$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a1, estimation.ueif.a0=crude.estimation$crude.cifj.sep.indirect.astar1.ate.ueif$ueif.a0,
                                                                                data.time=sort(unique(data.df$end.time)), estimation.time=sort(unique(data.df$end.time)))

relapse.separable.crude.ls$crude.rmtlj.sep.indirect.astar1.ate <- rmtlj_inference(npath=npath, train.ueif.a1=crude.train$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a1, train.ueif.a0=crude.train$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=crude.estimation$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a1, estimation.ueif.a0=crude.estimation$crude.rmtlj.sep.indirect.astar1.ate.ueif$ueif.a0,
                                                                                  data.time=sort(unique(data.df$end.time)), estimation.time=sort(unique(data.df$end.time)))
#
relapse.separable.crude.ls$crude.cifj.sep.direct.astar0.ate <- cifj_inference(npath=npath, train.ueif.a1=crude.train$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a1, train.ueif.a0=crude.train$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a0,
                                                                              estimation.ueif.a1=crude.estimation$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a1, estimation.ueif.a0=crude.estimation$crude.cifj.sep.direct.astar0.ate.ueif$ueif.a0,
                                                                              data.time=sort(unique(data.df$end.time)), estimation.time=sort(unique(data.df$end.time)))

relapse.separable.crude.ls$crude.rmtlj.sep.direct.astar0.ate <- rmtlj_inference(npath=npath, train.ueif.a1=crude.train$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a1, train.ueif.a0=crude.train$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a0,
                                                                                estimation.ueif.a1=crude.estimation$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a1, estimation.ueif.a0=crude.estimation$crude.rmtlj.sep.direct.astar0.ate.ueif$ueif.a0,
                                                                                data.time=sort(unique(data.df$end.time)), estimation.time=sort(unique(data.df$end.time)))
#
relapse.separable.crude.ls$crude.cifj.sep.indirect.astar0.ate <- cifj_inference(npath=npath, train.ueif.a1=crude.train$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a1, train.ueif.a0=crude.train$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a0,
                                                                                estimation.ueif.a1=crude.estimation$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a1, estimation.ueif.a0=crude.estimation$crude.cifj.sep.indirect.astar0.ate.ueif$ueif.a0,
                                                                                data.time=sort(unique(data.df$end.time)), estimation.time=sort(unique(data.df$end.time)))

relapse.separable.crude.ls$crude.rmtlj.sep.indirect.astar0.ate <- rmtlj_inference(npath=npath, train.ueif.a1=crude.train$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a1, train.ueif.a0=crude.train$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=crude.estimation$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a1, estimation.ueif.a0=crude.estimation$crude.rmtlj.sep.indirect.astar0.ate.ueif$ueif.a0,
                                                                                  data.time=sort(unique(data.df$end.time)), estimation.time=sort(unique(data.df$end.time)))

save(relapse.separable.crude.ls, file="relapse_separable_crude.RData")

