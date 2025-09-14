library(data.table)
library(dplyr)
library(table1)
library(xtable)

load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

####################################################################################
## cancer incidence curve
data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(cancer_eventdate)|(prescription_date<cancer_eventdate)) %>%
  filter(!(cancer==1 & is.na(cancer_eventdate))) %>%
  mutate(event1 = cancer) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(cancer_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(cancer, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & cancer==0), 2, cancer)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = cancer) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(cancer), cancer_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n, na.rm=TRUE))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(cancer), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(cancer), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n, na.rm=TRUE))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|cancer)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

data.df$a <- data.df$ref_sulf
table(data.df$ref_sulf)

## unweighted
# continuous
paste0(round(mean(data.df$prescription_age), digits=1), " (", round(sd(data.df$prescription_age), digits=1), ") & ",
       round(mean(data.df$prescription_age[data.df$ref_sulf==0]), digits=1), " (", round(sd(data.df$prescription_age[data.df$ref_sulf==0]), digits=1), ") & ",
       round(mean(data.df$prescription_age[data.df$ref_sulf==1]), digits=1), " (", round(sd(data.df$prescription_age[data.df$ref_sulf==1]), digits=1), ") & ",
       round(with(data.df, (mean(prescription_age[ref_sulf==1])-mean(prescription_age[ref_sulf==0]))/sqrt(var(prescription_age[ref_sulf==1])/sum(ref_sulf==1)+var(prescription_age[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(round(mean(data.df$prescription_year), digits=1), " (", round(sd(data.df$prescription_year), digits=1), ") & ",
       round(mean(data.df$prescription_year[data.df$ref_sulf==0]), digits=1), " (", round(sd(data.df$prescription_year[data.df$ref_sulf==0]), digits=1), ") & ",
       round(mean(data.df$prescription_year[data.df$ref_sulf==1]), digits=1), " (", round(sd(data.df$prescription_year[data.df$ref_sulf==1]), digits=1), ") & ",
       round(with(data.df, (mean(prescription_year[ref_sulf==1])-mean(prescription_year[ref_sulf==0]))/sqrt(var(prescription_year[ref_sulf==1])/sum(ref_sulf==1)+var(prescription_year[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(round(mean(data.df$bmi), digits=1), " (", round(sd(data.df$bmi), digits=1), ") & ",
       round(mean(data.df$bmi[data.df$ref_sulf==0]), digits=1), " (", round(sd(data.df$bmi[data.df$ref_sulf==0]), digits=1), ") & ",
       round(mean(data.df$bmi[data.df$ref_sulf==1]), digits=1), " (", round(sd(data.df$bmi[data.df$ref_sulf==1]), digits=1), ") & ",
       round(with(data.df, (mean(bmi[ref_sulf==1])-mean(bmi[ref_sulf==0]))/sqrt(var(bmi[ref_sulf==1])/sum(ref_sulf==1)+var(bmi[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(round(mean(data.df$hba1c_converge_1year_latest), digits=1), " (", round(sd(data.df$hba1c_converge_1year_latest), digits=1), ") & ",
       round(mean(data.df$hba1c_converge_1year_latest[data.df$ref_sulf==0]), digits=1), " (", round(sd(data.df$hba1c_converge_1year_latest[data.df$ref_sulf==0]), digits=1), ") & ",
       round(mean(data.df$hba1c_converge_1year_latest[data.df$ref_sulf==1]), digits=1), " (", round(sd(data.df$hba1c_converge_1year_latest[data.df$ref_sulf==1]), digits=1), ") & ",
       round(with(data.df, (mean(hba1c_converge_1year_latest[ref_sulf==1])-mean(hba1c_converge_1year_latest[ref_sulf==0]))/sqrt(var(hba1c_converge_1year_latest[ref_sulf==1])/sum(ref_sulf==1)+var(hba1c_converge_1year_latest[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

data.df$diabetes_after_prescription <- data.df$diabetes_after_prescription/365.25
paste0(round(mean(data.df$diabetes_after_prescription), digits=1), " (", round(sd(data.df$diabetes_after_prescription), digits=1), ") & ",
       round(mean(data.df$diabetes_after_prescription[data.df$ref_sulf==0]), digits=1), " (", round(sd(data.df$diabetes_after_prescription[data.df$ref_sulf==0]), digits=1), ") & ",
       round(mean(data.df$diabetes_after_prescription[data.df$ref_sulf==1]), digits=1), " (", round(sd(data.df$diabetes_after_prescription[data.df$ref_sulf==1]), digits=1), ") & ",
       round(with(data.df, (mean(diabetes_after_prescription[ref_sulf==1])-mean(diabetes_after_prescription[ref_sulf==0]))/sqrt(var(diabetes_after_prescription[ref_sulf==1])/sum(ref_sulf==1)+var(diabetes_after_prescription[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

# categorical
paste0(sum(data.df$gender), " & ", round(100*mean(data.df$gender), digits=1), "% & ",
       sum(data.df$gender[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$gender[data.df$ref_sulf==0]), digits=1), "% & ",
       sum(data.df$gender[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$gender[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean(gender[ref_sulf==1])-mean(gender[ref_sulf==0]))/sqrt(var(gender[ref_sulf==1])/sum(ref_sulf==1)+var(gender[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

# imd
paste0(sum(data.df$imd_imputed==1), " & ", round(100*mean(data.df$imd_imputed==1), digits=1), "% & ",
       sum((data.df$imd_imputed==1)[data.df$ref_sulf==0]), " & ", round(100*mean((data.df$imd_imputed==1)[data.df$ref_sulf==0]), digits=1), "% & ",
       sum((data.df$imd_imputed==1)[data.df$ref_sulf==1]), " & ", round(100*mean((data.df$imd_imputed==1)[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean((imd_imputed==1)[ref_sulf==1])-mean((imd_imputed==1)[ref_sulf==0]))/sqrt(var((imd_imputed==1)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==1)[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$imd_imputed==2), " & ", round(100*mean(data.df$imd_imputed==2), digits=1), "% & ",
       sum((data.df$imd_imputed==2)[data.df$ref_sulf==0]), " & ", round(100*mean((data.df$imd_imputed==2)[data.df$ref_sulf==0]), digits=1), "% & ",
       sum((data.df$imd_imputed==2)[data.df$ref_sulf==1]), " & ", round(100*mean((data.df$imd_imputed==2)[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean((imd_imputed==2)[ref_sulf==1])-mean((imd_imputed==2)[ref_sulf==0]))/sqrt(var((imd_imputed==2)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==2)[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$imd_imputed==3), " & ", round(100*mean(data.df$imd_imputed==3), digits=1), "% & ",
       sum((data.df$imd_imputed==3)[data.df$ref_sulf==0]), " & ", round(100*mean((data.df$imd_imputed==3)[data.df$ref_sulf==0]), digits=1), "% & ",
       sum((data.df$imd_imputed==3)[data.df$ref_sulf==1]), " & ", round(100*mean((data.df$imd_imputed==3)[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean((imd_imputed==3)[ref_sulf==1])-mean((imd_imputed==3)[ref_sulf==0]))/sqrt(var((imd_imputed==3)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==3)[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$imd_imputed==4), " & ", round(100*mean(data.df$imd_imputed==4), digits=1), "% & ",
       sum((data.df$imd_imputed==4)[data.df$ref_sulf==0]), " & ", round(100*mean((data.df$imd_imputed==4)[data.df$ref_sulf==0]), digits=1), "% & ",
       sum((data.df$imd_imputed==4)[data.df$ref_sulf==1]), " & ", round(100*mean((data.df$imd_imputed==4)[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean((imd_imputed==4)[ref_sulf==1])-mean((imd_imputed==4)[ref_sulf==0]))/sqrt(var((imd_imputed==4)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==4)[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$imd_imputed==5), " & ", round(100*mean(data.df$imd_imputed==5), digits=1), "% & ",
       sum((data.df$imd_imputed==5)[data.df$ref_sulf==0]), " & ", round(100*mean((data.df$imd_imputed==5)[data.df$ref_sulf==0]), digits=1), "% & ",
       sum((data.df$imd_imputed==5)[data.df$ref_sulf==1]), " & ", round(100*mean((data.df$imd_imputed==5)[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean((imd_imputed==5)[ref_sulf==1])-mean((imd_imputed==5)[ref_sulf==0]))/sqrt(var((imd_imputed==5)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==5)[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

# smoking
paste0(sum(data.df$smoke==0), " & ", round(100*mean(data.df$smoke==0), digits=1), "% & ",
       sum((data.df$smoke==0)[data.df$ref_sulf==0]), " & ", round(100*mean((data.df$smoke==0)[data.df$ref_sulf==0]), digits=1), "% & ",
       sum((data.df$smoke==0)[data.df$ref_sulf==1]), " & ", round(100*mean((data.df$smoke==0)[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean((smoke==0)[ref_sulf==1])-mean((smoke==0)[ref_sulf==0]))/sqrt(var((smoke==0)[ref_sulf==1])/sum(ref_sulf==1)+var((smoke==0)[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$smoke==1), " & ", round(100*mean(data.df$smoke==1), digits=1), "% & ",
       sum((data.df$smoke==1)[data.df$ref_sulf==0]), " & ", round(100*mean((data.df$smoke==1)[data.df$ref_sulf==0]), digits=1), "% & ",
       sum((data.df$smoke==1)[data.df$ref_sulf==1]), " & ", round(100*mean((data.df$smoke==1)[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean((smoke==1)[ref_sulf==1])-mean((smoke==1)[ref_sulf==0]))/sqrt(var((smoke==1)[ref_sulf==1])/sum(ref_sulf==1)+var((smoke==1)[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$smoke==2), " & ", round(100*mean(data.df$smoke==2), digits=1), "% & ",
       sum((data.df$smoke==2)[data.df$ref_sulf==0]), " & ", round(100*mean((data.df$smoke==2)[data.df$ref_sulf==0]), digits=1), "% & ",
       sum((data.df$smoke==2)[data.df$ref_sulf==1]), " & ", round(100*mean((data.df$smoke==2)[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean((smoke==2)[ref_sulf==1])-mean((smoke==2)[ref_sulf==0]))/sqrt(var((smoke==2)[ref_sulf==1])/sum(ref_sulf==1)+var((smoke==2)[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

# comorbidities
paste0(sum(data.df$hf), " & ", round(100*mean(data.df$hf), digits=1), "% & ",
       sum(data.df$hf[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$hf[data.df$ref_sulf==0]), digits=1), "% & ",
       sum(data.df$hf[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$hf[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean(hf[ref_sulf==1])-mean(hf[ref_sulf==0]))/sqrt(var(hf[ref_sulf==1])/sum(ref_sulf==1)+var(hf[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$chd_combined), " & ", round(100*mean(data.df$chd_combined), digits=1), "% & ",
       sum(data.df$chd_combined[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$chd_combined[data.df$ref_sulf==0]), digits=1), "% & ",
       sum(data.df$chd_combined[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$chd_combined[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean(chd_combined[ref_sulf==1])-mean(chd_combined[ref_sulf==0]))/sqrt(var(chd_combined[ref_sulf==1])/sum(ref_sulf==1)+var(chd_combined[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$af), " & ", round(100*mean(data.df$af), digits=1), "% & ",
       sum(data.df$af[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$af[data.df$ref_sulf==0]), digits=1), "% & ",
       sum(data.df$af[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$af[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean(af[ref_sulf==1])-mean(af[ref_sulf==0]))/sqrt(var(af[ref_sulf==1])/sum(ref_sulf==1)+var(af[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$stroke_combined), " & ", round(100*mean(data.df$stroke_combined), digits=1), "% & ",
       sum(data.df$stroke_combined[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$stroke_combined[data.df$ref_sulf==0]), digits=1), "% & ",
       sum(data.df$stroke_combined[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$stroke_combined[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean(stroke_combined[ref_sulf==1])-mean(stroke_combined[ref_sulf==0]))/sqrt(var(stroke_combined[ref_sulf==1])/sum(ref_sulf==1)+var(stroke_combined[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$hypertension_combined), " & ", round(100*mean(data.df$hypertension_combined), digits=1), "% & ",
       sum(data.df$hypertension_combined[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$hypertension_combined[data.df$ref_sulf==0]), digits=1), "% & ",
       sum(data.df$hypertension_combined[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$hypertension_combined[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean(hypertension_combined[ref_sulf==1])-mean(hypertension_combined[ref_sulf==0]))/sqrt(var(hypertension_combined[ref_sulf==1])/sum(ref_sulf==1)+var(hypertension_combined[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

paste0(sum(data.df$pvd), " & ", round(100*mean(data.df$pvd), digits=1), "% & ",
       sum(data.df$pvd[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$pvd[data.df$ref_sulf==0]), digits=1), "% & ",
       sum(data.df$pvd[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$pvd[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean(pvd[ref_sulf==1])-mean(pvd[ref_sulf==0]))/sqrt(var(pvd[ref_sulf==1])/sum(ref_sulf==1)+var(pvd[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

# paste0(sum(data.df$ckd), " & ", round(100*mean(data.df$ckd), digits=1), "% & ",
#        sum(data.df$ckd[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$ckd[data.df$ref_sulf==0]), digits=1), "% & ",
#        sum(data.df$ckd[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$ckd[data.df$ref_sulf==1]), digits=1), "%")

paste0(sum(data.df$copd), " & ", round(100*mean(data.df$copd), digits=1), "% & ",
       sum(data.df$copd[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$copd[data.df$ref_sulf==0]), digits=1), "% & ",
       sum(data.df$copd[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$copd[data.df$ref_sulf==1]), digits=1), "% & ",
       round(with(data.df, (mean(copd[ref_sulf==1])-mean(copd[ref_sulf==0]))/sqrt(var(copd[ref_sulf==1])/sum(ref_sulf==1)+var(copd[ref_sulf==0])/sum(ref_sulf==0))), digits=1))

## PS
a.sl.lib <- c("SL.mean","SL.glm","SL.nnet","SL.rpartPrune","SL.xgboost","SL.ranger","SL.step","SL.glmnet","SL.earth") # a.sl.lib <- c("SL.mean","SL.glm") # ,"SL.kernelKnn","SL.gam"
n <- nrow(data.df)
data.df$a <- data.df$ref_sulf

train.ind <- sample(1:n, ceiling(n/2), replace=FALSE) # sample training set
train.df <- data.df[train.ind,]
train.df <- train.df[order(train.df$patid),]
estimation.df <- data.df[setdiff(1:n, train.ind),]
estimation.df <- estimation.df[order(estimation.df$patid),] # estimation.df <- estimation.df[order(estimation.df$end.time,-estimation.df$event),]

## train
train.df$naive <- 1
train.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = estimation.df$a, X = estimation.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
train.df$iptw <- with(train.df, 1/(ps*a+(1-ps)*(1-a)))
train.df$ow <- with(train.df, ifelse(a,1-ps,ps))

## estimation
estimation.df$naive <- 1
estimation.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = train.df$a, X = train.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred) # estimation.df$ps <- predict(glm(a~x1+x2+x3+x4,data=estimation.df,family=binomial()),type="response")
estimation.df$iptw <- with(estimation.df, 1/(ps*a+(1-ps)*(1-a)))
estimation.df$ow <- with(estimation.df, ifelse(a,1-ps,ps))

data.df <- rbind(train.df, estimation.df)

## IPTW and OW
# continuous
with(data.df, paste0(round(sum(prescription_age*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(sum(prescription_age*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(prescription_age*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(prescription_age*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(prescription_age[ref_sulf==1])/sum(ref_sulf==1)+var(prescription_age[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(sum(prescription_age*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(sum(prescription_age*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(prescription_age*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(prescription_age*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(prescription_age[ref_sulf==1])/sum(ref_sulf==1)+var(prescription_age[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(sum(prescription_year*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(sum(prescription_year*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(prescription_year*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(prescription_year*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(prescription_year[ref_sulf==1])/sum(ref_sulf==1)+var(prescription_year[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(sum(prescription_year*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(sum(prescription_year*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(prescription_year*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(prescription_year*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(prescription_year[ref_sulf==1])/sum(ref_sulf==1)+var(prescription_year[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(sum(bmi*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(sum(bmi*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(bmi*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(bmi*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(bmi[ref_sulf==1])/sum(ref_sulf==1)+var(bmi[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(sum(bmi*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(sum(bmi*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(bmi*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(bmi*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(bmi[ref_sulf==1])/sum(ref_sulf==1)+var(bmi[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(sum(hba1c_converge_1year_latest*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(sum(hba1c_converge_1year_latest*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(hba1c_converge_1year_latest*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(hba1c_converge_1year_latest*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(hba1c_converge_1year_latest[ref_sulf==1])/sum(ref_sulf==1)+var(hba1c_converge_1year_latest[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(sum(hba1c_converge_1year_latest*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(sum(hba1c_converge_1year_latest*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(hba1c_converge_1year_latest*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(hba1c_converge_1year_latest*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(hba1c_converge_1year_latest[ref_sulf==1])/sum(ref_sulf==1)+var(hba1c_converge_1year_latest[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(sum(diabetes_after_prescription*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(sum(diabetes_after_prescription*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(diabetes_after_prescription*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(diabetes_after_prescription*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(diabetes_after_prescription[ref_sulf==1])/sum(ref_sulf==1)+var(diabetes_after_prescription[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(sum(diabetes_after_prescription*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(sum(diabetes_after_prescription*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(diabetes_after_prescription*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(diabetes_after_prescription*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(diabetes_after_prescription[ref_sulf==1])/sum(ref_sulf==1)+var(diabetes_after_prescription[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

# categorical
with(data.df, paste0(round(100*sum(gender*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum(gender*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(gender*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(gender*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(gender[ref_sulf==1])/sum(ref_sulf==1)+var(gender[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum(gender*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum(gender*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(gender*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(gender*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(gender[ref_sulf==1])/sum(ref_sulf==1)+var(gender[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

# imd
with(data.df, paste0(round(100*sum((imd_imputed==1)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum((imd_imputed==1)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum((imd_imputed==1)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum((imd_imputed==1)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var((imd_imputed==1)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==1)[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum((imd_imputed==1)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum((imd_imputed==1)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum((imd_imputed==1)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum((imd_imputed==1)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var((imd_imputed==1)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==1)[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum((imd_imputed==2)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum((imd_imputed==2)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum((imd_imputed==2)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum((imd_imputed==2)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var((imd_imputed==2)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==2)[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum((imd_imputed==2)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum((imd_imputed==2)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum((imd_imputed==2)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum((imd_imputed==2)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var((imd_imputed==2)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==2)[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum((imd_imputed==3)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum((imd_imputed==3)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum((imd_imputed==3)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum((imd_imputed==3)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var((imd_imputed==3)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==3)[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum((imd_imputed==3)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum((imd_imputed==3)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum((imd_imputed==3)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum((imd_imputed==3)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var((imd_imputed==3)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==3)[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum((imd_imputed==4)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum((imd_imputed==4)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum((imd_imputed==4)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum((imd_imputed==4)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var((imd_imputed==4)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==4)[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum((imd_imputed==4)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum((imd_imputed==4)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum((imd_imputed==4)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum((imd_imputed==4)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var((imd_imputed==4)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==4)[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum((imd_imputed==5)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum((imd_imputed==5)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum((imd_imputed==5)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum((imd_imputed==5)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var((imd_imputed==5)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==5)[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum((imd_imputed==5)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum((imd_imputed==5)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum((imd_imputed==5)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum((imd_imputed==5)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var((imd_imputed==5)[ref_sulf==1])/sum(ref_sulf==1)+var((imd_imputed==5)[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

# smoking
with(data.df, paste0(round(100*sum((smoke==0)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum((smoke==0)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum((smoke==0)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum((smoke==0)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var((smoke==0)[ref_sulf==1])/sum(ref_sulf==1)+var((smoke==0)[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum((smoke==0)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum((smoke==0)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum((smoke==0)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum((smoke==0)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var((smoke==0)[ref_sulf==1])/sum(ref_sulf==1)+var((smoke==0)[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum((smoke==1)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum((smoke==1)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum((smoke==1)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum((smoke==1)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var((smoke==1)[ref_sulf==1])/sum(ref_sulf==1)+var((smoke==1)[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum((smoke==1)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum((smoke==1)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum((smoke==1)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum((smoke==1)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var((smoke==1)[ref_sulf==1])/sum(ref_sulf==1)+var((smoke==1)[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum((smoke==2)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum((smoke==2)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum((smoke==2)*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum((smoke==2)*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var((smoke==2)[ref_sulf==1])/sum(ref_sulf==1)+var((smoke==2)[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum((smoke==2)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum((smoke==2)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum((smoke==2)*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum((smoke==2)*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var((smoke==2)[ref_sulf==1])/sum(ref_sulf==1)+var((smoke==2)[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

# comorbidities
with(data.df, paste0(round(100*sum(hf*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum(hf*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(hf*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(hf*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(hf[ref_sulf==1])/sum(ref_sulf==1)+var(hf[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum(hf*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum(hf*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(hf*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(hf*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(hf[ref_sulf==1])/sum(ref_sulf==1)+var(hf[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum(chd_combined*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum(chd_combined*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(chd_combined*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(chd_combined*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(chd_combined[ref_sulf==1])/sum(ref_sulf==1)+var(chd_combined[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum(chd_combined*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum(chd_combined*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(chd_combined*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(chd_combined*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(chd_combined[ref_sulf==1])/sum(ref_sulf==1)+var(chd_combined[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum(af*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum(af*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(af*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(af*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(af[ref_sulf==1])/sum(ref_sulf==1)+var(af[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum(af*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum(af*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(af*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(af*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(af[ref_sulf==1])/sum(ref_sulf==1)+var(af[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum(stroke_combined*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum(stroke_combined*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(stroke_combined*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(stroke_combined*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(stroke_combined[ref_sulf==1])/sum(ref_sulf==1)+var(stroke_combined[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum(stroke_combined*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum(stroke_combined*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(stroke_combined*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(stroke_combined*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(stroke_combined[ref_sulf==1])/sum(ref_sulf==1)+var(stroke_combined[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum(hypertension_combined*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum(hypertension_combined*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(hypertension_combined*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(hypertension_combined*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(hypertension_combined[ref_sulf==1])/sum(ref_sulf==1)+var(hypertension_combined[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum(hypertension_combined*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum(hypertension_combined*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(hypertension_combined*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(hypertension_combined*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(hypertension_combined[ref_sulf==1])/sum(ref_sulf==1)+var(hypertension_combined[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

with(data.df, paste0(round(100*sum(pvd*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum(pvd*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(pvd*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(pvd*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(pvd[ref_sulf==1])/sum(ref_sulf==1)+var(pvd[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum(pvd*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum(pvd*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(pvd*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(pvd*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(pvd[ref_sulf==1])/sum(ref_sulf==1)+var(pvd[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

# paste0(sum(data.df$ckd), " & ", round(100*mean(data.df$ckd), digits=1), "% & ",
#        sum(data.df$ckd[data.df$ref_sulf==0]), " & ", round(100*mean(data.df$ckd[data.df$ref_sulf==0]), digits=1), "% & ",
#        sum(data.df$ckd[data.df$ref_sulf==1]), " & ", round(100*mean(data.df$ckd[data.df$ref_sulf==1]), digits=1), "%")

with(data.df, paste0(round(100*sum(copd*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw), digits=1), " & ",
                     round(100*sum(copd*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw), digits=1), " & ",
                     round((sum(copd*(ref_sulf==1)*iptw)/sum((ref_sulf==1)*iptw)-sum(copd*(ref_sulf==0)*iptw)/sum((ref_sulf==0)*iptw))/sqrt(var(copd[ref_sulf==1])/sum(ref_sulf==1)+var(copd[ref_sulf==0])/sum(ref_sulf==0)), digits=1), " & ",
                     round(100*sum(copd*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow), digits=1), " & ",
                     round(100*sum(copd*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow), digits=1), " & ",
                     round((sum(copd*(ref_sulf==1)*ow)/sum((ref_sulf==1)*ow)-sum(copd*(ref_sulf==0)*ow)/sum((ref_sulf==0)*ow))/sqrt(var(copd[ref_sulf==1])/sum(ref_sulf==1)+var(copd[ref_sulf==0])/sum(ref_sulf==0)), digits=1)))

### cancer incidence

table(data.df$event3, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event3==0 & with(data.df, pmin(cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event3==0 & with(data.df, pmin(cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

### breast cancer incidence
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<breast_cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(breast_cancer_eventdate)|(prescription_date<breast_cancer_eventdate)) %>%
  filter(!(breast_cancer==1 & is.na(breast_cancer_eventdate))) %>%
  mutate(event1 = breast_cancer) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(breast_cancer_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(breast_cancer, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & breast_cancer==0), 2, breast_cancer)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = breast_cancer) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(breast_cancer), breast_cancer_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(breast_cancer), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(breast_cancer), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|breast_cancer)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

# gender==1 and removed from confounders
data.df <- data.df[data.df$gender==0, ]
confounders <- c("prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event3, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event3==0 & with(data.df, pmin(breast_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(breast_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event3==0 & with(data.df, pmin(breast_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(breast_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

### prostate cancer incidence
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<prostate_cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(prostate_cancer_eventdate)|(prescription_date<prostate_cancer_eventdate)) %>%
  filter(!(prostate_cancer==1 & is.na(prostate_cancer_eventdate))) %>%
  mutate(event1 = prostate_cancer) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(prostate_cancer_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(prostate_cancer, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & prostate_cancer==0), 2, prostate_cancer)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = prostate_cancer) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(prostate_cancer), prostate_cancer_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(prostate_cancer), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(prostate_cancer), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|prostate_cancer)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

data.df <- data.df[data.df$gender==1,]
confounders <- c("prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event3, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event3==0 & with(data.df, pmin(prostate_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(prostate_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event3==0 & with(data.df, pmin(prostate_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(prostate_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

### lung cancer incidence
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<lung_cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(lung_cancer_eventdate)|(prescription_date<lung_cancer_eventdate)) %>%
  filter(!(lung_cancer==1 & is.na(lung_cancer_eventdate))) %>%
  mutate(event1 = lung_cancer) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(lung_cancer_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(lung_cancer, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & lung_cancer==0), 2, lung_cancer)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = lung_cancer) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(lung_cancer), lung_cancer_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(lung_cancer), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(lung_cancer), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|lung_cancer)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event3, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event3==0 & with(data.df, pmin(lung_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(lung_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event3==0 & with(data.df, pmin(lung_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(lung_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

### bowel cancer incidence
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<bowel_cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(bowel_cancer_eventdate)|(prescription_date<bowel_cancer_eventdate)) %>%
  filter(!(bowel_cancer==1 & is.na(bowel_cancer_eventdate))) %>%
  mutate(event1 = bowel_cancer) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(bowel_cancer_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(bowel_cancer, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & bowel_cancer==0), 2, bowel_cancer)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = bowel_cancer) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(bowel_cancer), bowel_cancer_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(bowel_cancer), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(bowel_cancer), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|bowel_cancer)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event3, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event3==0 & with(data.df, pmin(bowel_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(bowel_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event3==0 & with(data.df, pmin(bowel_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(bowel_cancer_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

### leukemia incidence
load("antidiabetes_aje.RData")
load("../AJE/subtype_extract.RData")
data.df <- merge(data.df, leukaemia.df, by="patid", all.x=TRUE)
data.df$leukaemia <- ifelse(is.na(data.df$leukaemia_eventdate), 0, 1)

data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<leukaemia_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(leukaemia_eventdate)|(prescription_date<leukaemia_eventdate)) %>%
  filter(!(leukaemia==1 & is.na(leukaemia_eventdate))) %>%
  mutate(event1 = leukaemia) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(leukaemia_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(leukaemia, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & leukaemia==0), 2, leukaemia)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = leukaemia) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(leukaemia), leukaemia_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(leukaemia), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(leukaemia), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|leukaemia)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event3, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event3==0 & with(data.df, pmin(leukaemia_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(leukaemia_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event3==0 & with(data.df, pmin(leukaemia_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(leukaemia_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

## lymphoma incidence
load("antidiabetes_aje.RData")
load("../AJE/subtype_extract.RData")
data.df <- merge(data.df, lymphoma.df, by="patid", all.x=TRUE)
data.df$lymphoma <- ifelse(is.na(data.df$lymphoma_eventdate), 0, 1)

data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<lymphoma_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(lymphoma_eventdate)|(prescription_date<lymphoma_eventdate)) %>%
  filter(!(lymphoma==1 & is.na(lymphoma_eventdate))) %>%
  mutate(event1 = lymphoma) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(lymphoma_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(lymphoma, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & lymphoma==0), 2, lymphoma)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = lymphoma) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(lymphoma), lymphoma_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(lymphoma), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(lymphoma), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|lymphoma)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event3, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event3==0 & with(data.df, pmin(lymphoma_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(lymphoma_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event3==0 & with(data.df, pmin(lymphoma_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(lymphoma_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

## myeloma incidence
load("antidiabetes_aje.RData")
load("../AJE/subtype_extract.RData")
data.df <- merge(data.df, myeloma.df, by="patid", all.x=TRUE)
data.df$myeloma <- ifelse(is.na(data.df$myeloma_eventdate), 0, 1)

data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<myeloma_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(myeloma_eventdate)|(prescription_date<myeloma_eventdate)) %>%
  filter(!(myeloma==1 & is.na(myeloma_eventdate))) %>%
  mutate(event1 = myeloma) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(myeloma_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(myeloma, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & myeloma==0), 2, myeloma)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = myeloma) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(myeloma), myeloma_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(myeloma), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(myeloma), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|myeloma)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event3, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event3==0 & with(data.df, pmin(myeloma_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(myeloma_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event3==0 & with(data.df, pmin(myeloma_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event3==0 & with(data.df, pmin(myeloma_eventdate, deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

### natural mortality
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>%
  mutate(event1 = death) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event1, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

### cancer mortality
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>%
  filter(prescription_date>cancer_eventdate) %>%
  mutate(event1 = death) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event1, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)

## breast cancer mortality
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<breast_cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>%
  filter(prescription_date>breast_cancer_eventdate) %>%
  mutate(event1 = death) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  data.frame()

# gender==0 and removed from confounders
data.df <- data.df[data.df$gender==0,]
confounders <- c("prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]

table(data.df$event1, data.df$ref_sulf)
table(data.df$ref_sulf)

sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))==data.df$lcd_n & data.df$ref_sulf==1)

sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==0)
sum(data.df$event1==0 & with(data.df, pmin(deathdate_n, lcd_n, tod_n, na.rm=TRUE))!=data.df$lcd_n & data.df$ref_sulf==1)









