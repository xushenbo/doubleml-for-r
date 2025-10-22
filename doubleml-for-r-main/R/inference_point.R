# Helper Functions for Survival Analysis Inference
# These functions support both survival analysis and competing risks methods

#' Last Observation Carried Forward
#' 
#' Simple implementation of na.locf (Last Observation Carried Forward)
#' Similar to zoo::na.locf functionality
#' 
#' @param x Vector or matrix to fill missing values
#' @return Object with missing values filled using last observation
na.locf <- function(x) {
  if (is.vector(x)) {
    for (i in 2:length(x)) {
      if (is.na(x[i])) {
        x[i] <- x[i-1]
      }
    }
  } else if (is.matrix(x)) {
    for (j in 1:ncol(x)) {
      for (i in 2:nrow(x)) {
        if (is.na(x[i, j])) {
          x[i, j] <- x[i-1, j]
        }
      }
    }
  }
  return(x)
}

#' Time-based Approximation for Survival Functions
#' 
#' Performs linear interpolation for time-dependent survival functions
#' Similar to approx() but with specific rules for survival data
#' 
#' @param object Data frame with time column and value columns
#' @param xout Output time points for interpolation
#' @param rule Integer specifying how to handle out-of-bounds values
#' @return Data frame with interpolated values
approxTime <- function(object, xout, rule = 2) {
  time_col <- object[, 1]  # First column is time
  result <- data.frame(t = xout)
  
  # For each non-time column, perform interpolation
  for (j in 2:ncol(object)) {
    col_values <- object[, j]
    # Use approx for linear interpolation
    interpolated <- approx(x = time_col, y = col_values, xout = xout, rule = rule)
    result[, j] <- interpolated$y
  }
  
  # Preserve column names
  names(result) <- names(object)
  return(result)
}

#' Extract Average Treatment Effect for Survival Probability at Tau
#'
#' Extracts point estimates, standard errors, and treatment effects for survival
#' probability at a specific time point tau using cumulative minimum for
#' monotonic decreasing survival probabilities.
#'
#' @param train.ueif.a1 UEIF matrix for treatment a=1 from training set
#' @param train.ueif.a0 UEIF matrix for treatment a=0 from training set
#' @param estimation.ueif.a1 UEIF matrix for treatment a=1 from estimation set
#' @param estimation.ueif.a0 UEIF matrix for treatment a=0 from estimation set
#' @param train.time Time points from training set
#' @param estimation.time Time points from estimation set
#' @param tau Time point at which to extract the treatment effect
#' @return Data frame with point estimates and standard errors at tau
survival_prob_extract_average <- function(train.ueif.a1, train.ueif.a0, estimation.ueif.a1, estimation.ueif.a0, train.time, estimation.time, tau){
  ## point estimate
  # a=0
  ueif.a0 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"), all=TRUE, by="t")
  ueif.a0 <- approxTime(ueif.a0, xout=ueif.a0$t)
  if(sum(is.na(ueif.a0))>0)
  {ueif.a0 <- approxTime(rbind(0, ueif.a0), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a0 <- t(ueif.a0[,-1])

  # Calculate point estimates with cummin for monotonic decreasing survival probabilities
  estimate.a0.raw <- colMeans(ueif.a0, na.rm = TRUE)
  estimate.a0 <- cummin(estimate.a0.raw)
  ueif.df <- data.frame(t=sort(unique(c(train.time, estimation.time)))[1:length(estimate.a0)], estimate.a0=estimate.a0)

  # a=1
  ueif.a1 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t") , all=TRUE, by="t")
  ueif.a1 <- approxTime(ueif.a1, xout=ueif.a1$t)
  if(sum(is.na(ueif.a1))>0)
  {ueif.a1 <- approxTime(rbind(0, ueif.a1), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a1 <- t(ueif.a1[,-1])

  # Calculate point estimates with cummin for monotonic decreasing survival probabilities
  estimate.a1.raw <- colMeans(ueif.a1, na.rm = TRUE)
  estimate.a1 <- cummin(estimate.a1.raw)
  ueif.df$estimate.a1 <- estimate.a1[1:nrow(ueif.df)]

  ueif.df <- approxTime(ueif.df, xout=ueif.df$t)
  ueif.df$estimate.diff <- ueif.df$estimate.a1-ueif.df$estimate.a0

  ## standard error - using cummax for monotonic increasing variance over time
  ueif.df$estimate.a0.se <- cummax(sqrt(colMeans((ueif.a0-matrix(ueif.df$estimate.a0, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(nrow(ueif.a0)))
  ueif.df$estimate.a1.se <- cummax(sqrt(colMeans((ueif.a1-matrix(ueif.df$estimate.a1, nrow=nrow(ueif.a1), ncol=ncol(ueif.a1), byrow=TRUE))^2))/sqrt(nrow(ueif.a1)))
  ueif.df$estimate.diff.se <- cummax(sqrt(colMeans(((ueif.a1-ueif.a0)-matrix(ueif.df$estimate.diff, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(nrow(ueif.a0)))

  # Use <= for right-continuous step function (cadlag)
  tau.ind <- max(which(ueif.df$t <= tau))
  ueif.df <- ueif.df[tau.ind,]
  return(ueif.df[,-1])
}

#' Extract Average Treatment Effect for Cumulative Incidence at Tau
#'
#' Extracts point estimates, standard errors, and treatment effects for cumulative
#' incidence function (CIF) at a specific time point tau using cumulative maximum
#' for monotonic increasing cumulative incidence probabilities.
#'
#' @param train.ueif.a1 UEIF matrix for treatment a=1 from training set
#' @param train.ueif.a0 UEIF matrix for treatment a=0 from training set
#' @param estimation.ueif.a1 UEIF matrix for treatment a=1 from estimation set
#' @param estimation.ueif.a0 UEIF matrix for treatment a=0 from estimation set
#' @param train.time Time points from training set
#' @param estimation.time Time points from estimation set
#' @param tau Time point at which to extract the treatment effect
#' @return Data frame with point estimates and standard errors at tau
cifj_extract_average <- function(train.ueif.a1, train.ueif.a0, estimation.ueif.a1, estimation.ueif.a0, train.time, estimation.time, tau){
  ## point estimate
  # a=0
  ueif.a0 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"), all=TRUE, by="t")
  ueif.a0 <- approxTime(ueif.a0, xout=ueif.a0$t)
  if(sum(is.na(ueif.a0))>0)
  {ueif.a0 <- approxTime(rbind(0, ueif.a0), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a0 <- t(ueif.a0[,-1])

  # Calculate point estimates with cummax for monotonic increasing CIF (no gcmlcm)
  estimate.a0.raw <- colMeans(ueif.a0, na.rm = TRUE)
  estimate.a0 <- cummax(estimate.a0.raw)
  ueif.df <- data.frame(t=sort(unique(c(train.time, estimation.time)))[1:length(estimate.a0)], estimate.a0=estimate.a0)

  # a=1
  ueif.a1 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t") , all=TRUE, by="t")
  ueif.a1 <- approxTime(ueif.a1, xout=ueif.a1$t)
  if(sum(is.na(ueif.a1))>0)
  {ueif.a1 <- approxTime(rbind(0, ueif.a1), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a1 <- t(ueif.a1[,-1])

  # Calculate point estimates with cummax for monotonic increasing CIF (no gcmlcm)
  estimate.a1.raw <- colMeans(ueif.a1, na.rm = TRUE)
  estimate.a1 <- cummax(estimate.a1.raw)
  ueif.df$estimate.a1 <- estimate.a1[1:nrow(ueif.df)]

  ueif.df <- approxTime(ueif.df, xout=ueif.df$t)
  ueif.df$estimate.diff <- ueif.df$estimate.a1-ueif.df$estimate.a0

  ## standard error - using cummax for monotonic increasing variance over time
  ueif.df$estimate.a0.se <- cummax(sqrt(colMeans((ueif.a0-matrix(ueif.df$estimate.a0, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(nrow(ueif.a0)))
  ueif.df$estimate.a1.se <- cummax(sqrt(colMeans((ueif.a1-matrix(ueif.df$estimate.a1, nrow=nrow(ueif.a1), ncol=ncol(ueif.a1), byrow=TRUE))^2))/sqrt(nrow(ueif.a1)))
  ueif.df$estimate.diff.se <- cummax(sqrt(colMeans(((ueif.a1-ueif.a0)-matrix(ueif.df$estimate.diff, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(nrow(ueif.a0)))

  # Use < for left-continuous behavior, no interpolation
  tau.ind <- max(which(ueif.df$t < tau))
  ueif.df <- ueif.df[tau.ind,]
  return(ueif.df[,-1])
}

#' Extract Average Treatment Effect for RMST at Tau
#'
#' Extracts point estimates, standard errors, and treatment effects for restricted
#' mean survival time (RMST) at a specific time point tau using least concave majorant
#' for monotonic increasing RMST values.
#'
#' @param train.ueif.a1 UEIF matrix for treatment a=1 from training set
#' @param train.ueif.a0 UEIF matrix for treatment a=0 from training set
#' @param estimation.ueif.a1 UEIF matrix for treatment a=1 from estimation set
#' @param estimation.ueif.a0 UEIF matrix for treatment a=0 from estimation set
#' @param train.time Time points from training set
#' @param estimation.time Time points from estimation set
#' @param tau Time point at which to extract the treatment effect
#' @return Data frame with point estimates and standard errors at tau
rmst_extract_average <- function(train.ueif.a1, train.ueif.a0, estimation.ueif.a1, estimation.ueif.a0, train.time, estimation.time, tau){
  ## point estimate
  # a=0
  ueif.a0 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"), all=TRUE, by="t")
  ueif.a0 <- approxTime(ueif.a0, xout=ueif.a0$t)
  if(sum(is.na(ueif.a0))>0)
  {ueif.a0 <- approxTime(rbind(0, ueif.a0), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a0 <- t(ueif.a0[,-1])

  estimate.a0.lcm <- gcmlcm(x=sort(unique(c(train.time, estimation.time)))[1:max(which(!is.na(colMeans(ueif.a0))))], y=colMeans(ueif.a0)[1:max(which(!is.na(colMeans(ueif.a0))))], type="lcm")[c("x.knots","y.knots")]
  names(estimate.a0.lcm) <- c("t", "estimate.a0")
  ueif.df <- merge(data.frame(t=sort(unique(c(train.time, estimation.time)))), estimate.a0.lcm, all=TRUE, by="t")

  # a=1
  ueif.a1 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t") , all=TRUE, by="t")
  ueif.a1 <- approxTime(ueif.a1, xout=ueif.a1$t)
  if(sum(is.na(ueif.a1))>0)
  {ueif.a1 <- approxTime(rbind(0, ueif.a1), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a1 <- t(ueif.a1[,-1])
  estimate.a1.lcm <- gcmlcm(x=sort(unique(c(train.time, estimation.time)))[1:max(which(!is.na(colMeans(ueif.a1))))], y=colMeans(ueif.a1)[1:max(which(!is.na(colMeans(ueif.a1))))], type="lcm")[c("x.knots","y.knots")]
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

#' Extract Average Treatment Effect for RMTL at Tau
#'
#' Extracts point estimates, standard errors, and treatment effects for restricted
#' mean time lost (RMTL) at a specific time point tau using greatest convex minorant
#' for monotonic increasing RMTL values.
#'
#' @param train.ueif.a1 UEIF matrix for treatment a=1 from training set
#' @param train.ueif.a0 UEIF matrix for treatment a=0 from training set
#' @param estimation.ueif.a1 UEIF matrix for treatment a=1 from estimation set
#' @param estimation.ueif.a0 UEIF matrix for treatment a=0 from estimation set
#' @param train.time Time points from training set
#' @param estimation.time Time points from estimation set
#' @param tau Time point at which to extract the treatment effect
#' @return Data frame with point estimates and standard errors at tau
rmtlj_extract_average <- function(train.ueif.a1, train.ueif.a0, estimation.ueif.a1, estimation.ueif.a0, train.time, estimation.time, tau){
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