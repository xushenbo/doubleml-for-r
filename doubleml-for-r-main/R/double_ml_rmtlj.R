#' @title Restricted Mean Time Lost for Competing Risks (RMTLJ)
#'
#' @description
#' Computes efficient influence function for restricted mean time lost to cause j
#' before time tau in competing risks settings.
#'
#' @param id (`numeric`) \\cr
#' Individual identifiers for observations.
#'
#' @param a (`numeric`) \\cr
#' Binary treatment indicator (0 or 1).
#'
#' @param time (`numeric`) \\cr
#' Observed time (minimum of event time and censoring time).
#'
#' @param event (`numeric`) \\cr
#' Event indicator (0 = censored, j = cause j event).
#'
#' @param bw (`numeric`) \\cr
#' Balancing weights (inverse propensity weights).
#'
#' @param tilt (`numeric`) \\cr
#' Tilting function for stabilization.
#'
#' @param G.a0, G.a1 (`matrix`) \\cr
#' Censoring survival function G(t|X,A=0) and G(t|X,A=1).
#'
#' @param S.a0, S.a1 (`matrix`) \\cr
#' Overall event-free survival function S(t|X,A=0) and S(t|X,A=1).
#'
#' @param Sj.a0, Sj.a1 (`matrix`) \\cr
#' Cause-specific survival function for cause j: Sj(t|X,A=0) and Sj(t|X,A=1).
#'
#' @param freq.time (`numeric(1)`) \\cr
#' Time grid frequency. If NULL, uses unique event times.
#'
#' @param admin.cens (`numeric(1)`) \\cr
#' Administrative censoring time (maximum follow-up).
#'
#' @return List with elements:
#' \describe{
#'   \item{ueif.a1}{Efficient influence function for A=1}
#'   \item{ueif.a0}{Efficient influence function for A=0}
#' }
#'
#' @details
#' This function implements the censoring martingale based efficient influence
#' function for restricted mean time lost (RMTL) to cause j in competing risks.
#' The RMTL measures the expected time lost to a specific cause before time tau.
#'
#' The efficient influence function has 5 terms:
#' \enumerate{
#'   \item Inverse probability weighted observed RMTL
#'   \item Predicted RMTL under treatment level
#'   \item Censoring correction for cumulative incidence
#'   \item Survival correction for RMTL
#'   \item Combined censoring and survival correction
#' }
#'
#' @keywords internal
double_ml_rmtlj <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))

  # Winsorization of survival functions
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))

  # Compute hazards
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))

  # Cumulative incidence and RMTL
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))

  # Indicators
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  RMTLj.obs <- do.call(cbind, lapply(1:ns, function(u)  (s[u]-pmin(time, s[u]))*(event == cause)))

  # 5-term UEIF for A=0
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.obs/matrix(rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0*(dNct - Yt*G.dHazard.a0)/(G.a0*S.a0), 1, cumsum))-
                                                              matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0*(dNct - Yt*G.dHazard.a0)/(G.a0*S.a0), 1, cumsum)))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0*(t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))-1)
  term5.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMTLj.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  ueif.a0 <- 1/mean(tilt)*(term1.a0+term2.a0+term3.a0+term4.a0-term5.a0)

  # 5-term UEIF for A=1 (symmetric)
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.obs/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1*(dNct - Yt*G.dHazard.a1)/(G.a1*S.a1), 1, cumsum))-
                                                              matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1*(dNct - Yt*G.dHazard.a1)/(G.a1*S.a1), 1, cumsum)))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term5.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMTLj.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  ueif.a1 <- 1/mean(tilt)*(term1.a1+term2.a1+term3.a1+term4.a1-term5.a1)

  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

#' @title Separable Direct Effect for RMTLJ with A*=1
#'
#' @description
#' Computes efficient influence function for separable direct effect of RMTL
#' with reference treatment A*=1 in competing risks mediation analysis.
#'
#' @param id (`numeric`) \\cr
#' Individual identifiers for observations.
#'
#' @param a (`numeric`) \\cr
#' Binary treatment indicator (0 or 1).
#'
#' @param time (`numeric`) \\cr
#' Observed time (minimum of event time and censoring time).
#'
#' @param event (`numeric`) \\cr
#' Event indicator (0 = censored, j = cause j event).
#'
#' @param bw (`numeric`) \\cr
#' Balancing weights (inverse propensity weights).
#'
#' @param tilt (`numeric`) \\cr
#' Tilting function for stabilization.
#'
#' @param G.a0, G.a1 (`matrix`) \\cr
#' Censoring survival function G(t|X,A=0) and G(t|X,A=1).
#'
#' @param S.a0, S.a1 (`matrix`) \\cr
#' Overall event-free survival function S(t|X,A=0) and S(t|X,A=1).
#'
#' @param Sj.a0, Sj.a1 (`matrix`) \\cr
#' Cause-specific survival for cause j: Sj(t|X,A=0) and Sj(t|X,A=1).
#'
#' @param Sjbar.a0, Sjbar.a1 (`matrix`) \\cr
#' Complementary cause survival (all causes except j): Sjbar(t|X,A=0) and Sjbar(t|X,A=1).
#'
#' @param freq.time (`numeric(1)`) \\cr
#' Time grid frequency. If NULL, uses unique event times.
#'
#' @param admin.cens (`numeric(1)`) \\cr
#' Administrative censoring time (maximum follow-up).
#'
#' @return List with elements:
#' \describe{
#'   \item{ueif.a1}{Efficient influence function for A=1 under separable direct effect}
#'   \item{ueif.a0}{Efficient influence function for A=0 under separable direct effect}
#' }
#'
#' @details
#' This function implements the separable direct effect decomposition where the
#' treatment affects cause j directly while holding the effect on complementary
#' causes at reference level A*=1. The counterfactual RMTL combines Sj under A=0
#' with Sjbar under A=1.
#'
#' The efficient influence function has 6 terms accounting for:
#' \enumerate{
#'   \item Inverse probability weighted observed RMTL
#'   \item Predicted counterfactual RMTL
#'   \item Censoring correction for cumulative incidence
#'   \item Survival correction for RMTL
#'   \item Combined censoring and survival correction
#'   \item Complementary cause correction term
#' }
#'
#' @keywords internal
double_ml_rmtlj_sep_direct_astar1 <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))

  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Sjbar.a0 <- t(na.locf(t(ifelse(Sjbar.a0 < 1e-3, 1e-3, Sjbar.a0))))
  Sjbar.a1 <- t(na.locf(t(ifelse(Sjbar.a1 < 1e-3, 1e-3, Sjbar.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fjbar.dHazard.a0 <- t(apply(cbind(0, -log(Sjbar.a0)), 1, diff))
  Fjbar.dHazard.a1 <- t(apply(cbind(0, -log(Sjbar.a1)), 1, diff))

  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  Fj.a0.a1 <- t(apply(cbind(1, winsorize_values(Sj.a0*Sjbar.a1, lower = 1e-3, upper = 1)[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a1!

  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a0.a1 <- t(apply(Fj.a0.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))

  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))

  # a=0
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0.a1
  term2.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a1*(dNjt-Yt*Fj.dHazard.a0)/(Sjbar.a0*G.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a1*(dNjt-Yt*Fj.dHazard.a0)/(Sjbar.a0*G.a0), 1, cumsum)))
  term3.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a0*t(apply(Sjbar.a1*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum))
                                                             -t(apply(RMTLj.a0*Sjbar.a1*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum)))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a1*Fj.a0*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a1*Fj.a0*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum)))
  term5.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                            -t(apply(RMTLj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  term6.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                             -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  ueif.a0 <- 1/mean(bw*(a==0))*(term2.a0+term3.a0+term4.a0+term5.a0+term6.a0)+1/mean(tilt)*term1.a0

  # a=1
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1
  term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a1*(dNjt-Yt*Fj.dHazard.a1)/(Sjbar.a1*G.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a1*(dNjt-Yt*Fj.dHazard.a1)/(Sjbar.a1*G.a1), 1, cumsum)))
  term3.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a1*t(apply(Sjbar.a1*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum))
                                                             -t(apply(RMTLj.a1*Sjbar.a1*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum)))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a1*Fj.a1*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a1*Fj.a1*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum)))
  term5.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                            -t(apply(RMTLj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  term6.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                             -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  ueif.a1 <- 1/mean(bw*(a==1))*(term2.a1+term3.a1+term4.a1+term5.a1+term6.a1)+1/mean(tilt)*term1.a1

  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

#' @title Separable Indirect Effect for RMTLJ with A*=1
#'
#' @description
#' Computes efficient influence function for separable indirect effect of RMTL
#' with reference treatment A*=1 in competing risks mediation analysis.
#'
#' @param id (`numeric`) \\cr
#' Individual identifiers for observations.
#'
#' @param a (`numeric`) \\cr
#' Binary treatment indicator (0 or 1).
#'
#' @param time (`numeric`) \\cr
#' Observed time (minimum of event time and censoring time).
#'
#' @param event (`numeric`) \\cr
#' Event indicator (0 = censored, j = cause j event).
#'
#' @param bw (`numeric`) \\cr
#' Balancing weights (inverse propensity weights).
#'
#' @param tilt (`numeric`) \\cr
#' Tilting function for stabilization.
#'
#' @param G.a0, G.a1 (`matrix`) \\cr
#' Censoring survival function G(t|X,A=0) and G(t|X,A=1).
#'
#' @param S.a0, S.a1 (`matrix`) \\cr
#' Overall event-free survival function S(t|X,A=0) and S(t|X,A=1).
#'
#' @param Sj.a0, Sj.a1 (`matrix`) \\cr
#' Cause-specific survival for cause j: Sj(t|X,A=0) and Sj(t|X,A=1).
#'
#' @param Sjbar.a0, Sjbar.a1 (`matrix`) \\cr
#' Complementary cause survival (all causes except j): Sjbar(t|X,A=0) and Sjbar(t|X,A=1).
#'
#' @param freq.time (`numeric(1)`) \\cr
#' Time grid frequency. If NULL, uses unique event times.
#'
#' @param admin.cens (`numeric(1)`) \\cr
#' Administrative censoring time (maximum follow-up).
#'
#' @return List with elements:
#' \describe{
#'   \item{ueif.a1}{Efficient influence function for A=1 under separable indirect effect}
#'   \item{ueif.a0}{Efficient influence function for A=0 under separable indirect effect}
#' }
#'
#' @details
#' This function implements the separable indirect effect decomposition where the
#' treatment affects complementary causes while holding the effect on cause j at
#' reference level A*=1. The counterfactual RMTL combines Sj under A=1 with Sjbar
#' under the treatment level.
#'
#' The efficient influence function has 6 terms similar to the direct effect but
#' with different complementary cause weighting.
#'
#' @keywords internal
double_ml_rmtlj_sep_indirect_astar1 <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))

  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Sjbar.a0 <- t(na.locf(t(ifelse(Sjbar.a0 < 1e-3, 1e-3, Sjbar.a0))))
  Sjbar.a1 <- t(na.locf(t(ifelse(Sjbar.a1 < 1e-3, 1e-3, Sjbar.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fjbar.dHazard.a0 <- t(apply(cbind(0, -log(Sjbar.a0)), 1, diff))
  Fjbar.dHazard.a1 <- t(apply(cbind(0, -log(Sjbar.a1)), 1, diff))

  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  Fj.a0.a1 <- t(apply(cbind(1, winsorize_values(Sj.a0*Sjbar.a1, lower = 1e-3, upper = 1)[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a1!

  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a0.a1 <- t(apply(Fj.a0.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))

  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))

  # a=0
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0
  term2.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply((dNjt-Yt*Fj.dHazard.a0)/(G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNjt-Yt*Fj.dHazard.a0)/(G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum)))
  term3.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a0*t(apply((dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum))
                                                             -t(apply(RMTLj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum)))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a1/Sjbar.a0), 1, cumsum)))
  term5.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                             -t(apply(RMTLj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  term6.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  ueif.a0 <- 1/mean(bw*(a==0))*(term2.a0+term3.a0+term4.a0+term5.a0+term6.a0)+1/mean(tilt)*term1.a0

  # a=1
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0.a1
  term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply((dNjt-Yt*Fj.dHazard.a1)/(G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNjt-Yt*Fj.dHazard.a1)/(G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum)))
  term3.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a1*t(apply((dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum))
                                                             -t(apply(RMTLj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum)))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a1/Sjbar.a1), 1, cumsum)))
  term5.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                             -t(apply(RMTLj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  term6.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  ueif.a1 <- 1/mean(bw*(a==1))*(term2.a1+term3.a1+term4.a1+term5.a1+term6.a1)+1/mean(tilt)*term1.a1

  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

#' @title Separable Direct Effect for RMTLJ with A*=0
#'
#' @description
#' Computes efficient influence function for separable direct effect of RMTL
#' with reference treatment A*=0 in competing risks mediation analysis.
#'
#' @param id (`numeric`) \\cr
#' Individual identifiers for observations.
#'
#' @param a (`numeric`) \\cr
#' Binary treatment indicator (0 or 1).
#'
#' @param time (`numeric`) \\cr
#' Observed time (minimum of event time and censoring time).
#'
#' @param event (`numeric`) \\cr
#' Event indicator (0 = censored, j = cause j event).
#'
#' @param bw (`numeric`) \\cr
#' Balancing weights (inverse propensity weights).
#'
#' @param tilt (`numeric`) \\cr
#' Tilting function for stabilization.
#'
#' @param G.a0, G.a1 (`matrix`) \\cr
#' Censoring survival function G(t|X,A=0) and G(t|X,A=1).
#'
#' @param S.a0, S.a1 (`matrix`) \\cr
#' Overall event-free survival function S(t|X,A=0) and S(t|X,A=1).
#'
#' @param Sj.a0, Sj.a1 (`matrix`) \\cr
#' Cause-specific survival for cause j: Sj(t|X,A=0) and Sj(t|X,A=1).
#'
#' @param Sjbar.a0, Sjbar.a1 (`matrix`) \\cr
#' Complementary cause survival (all causes except j): Sjbar(t|X,A=0) and Sjbar(t|X,A=1).
#'
#' @param freq.time (`numeric(1)`) \\cr
#' Time grid frequency. If NULL, uses unique event times.
#'
#' @param admin.cens (`numeric(1)`) \\cr
#' Administrative censoring time (maximum follow-up).
#'
#' @return List with elements:
#' \describe{
#'   \item{ueif.a1}{Efficient influence function for A=1 under separable direct effect}
#'   \item{ueif.a0}{Efficient influence function for A=0 under separable direct effect}
#' }
#'
#' @details
#' This function implements the separable direct effect decomposition where the
#' treatment affects cause j directly while holding the effect on complementary
#' causes at reference level A*=0. The counterfactual RMTL combines Sj under A=1
#' with Sjbar under A=0.
#'
#' The efficient influence function has 6 terms, symmetric to the A*=1 version but
#' with reference level A*=0.
#'
#' @keywords internal
double_ml_rmtlj_sep_direct_astar0 <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))

  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Sjbar.a0 <- t(na.locf(t(ifelse(Sjbar.a0 < 1e-3, 1e-3, Sjbar.a0))))
  Sjbar.a1 <- t(na.locf(t(ifelse(Sjbar.a1 < 1e-3, 1e-3, Sjbar.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fjbar.dHazard.a0 <- t(apply(cbind(0, -log(Sjbar.a0)), 1, diff))
  Fjbar.dHazard.a1 <- t(apply(cbind(0, -log(Sjbar.a1)), 1, diff))

  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  Fj.a1.a0 <- t(apply(cbind(1, winsorize_values(Sj.a1*Sjbar.a0, lower = 1e-3, upper = 1)[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum)) # not Fj.dHazard.a1!

  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1.a0 <- t(apply(Fj.a1.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))

  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))

  # a=0
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0
  term2.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a0*(dNjt-Yt*Fj.dHazard.a0)/(Sjbar.a0*G.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a0*(dNjt-Yt*Fj.dHazard.a0)/(Sjbar.a0*G.a0), 1, cumsum)))
  term3.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a0*t(apply(Sjbar.a0*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum))
                                                             -t(apply(RMTLj.a0*Sjbar.a0*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum)))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a0*Fj.a0*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a0*Fj.a0*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum)))
  term5.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                            -t(apply(RMTLj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  term6.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                             -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  ueif.a0 <- 1/mean(bw*(a==0))*(term2.a0+term3.a0+term4.a0+term5.a0+term6.a0)+1/mean(tilt)*term1.a0

  # a=1
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1.a0
  term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a0*(dNjt-Yt*Fj.dHazard.a1)/(Sjbar.a1*G.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a0*(dNjt-Yt*Fj.dHazard.a1)/(Sjbar.a1*G.a1), 1, cumsum)))
  term3.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a1*t(apply(Sjbar.a0*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum))
                                                             -t(apply(RMTLj.a1*Sjbar.a0*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum)))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Sjbar.a0*Fj.a1*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Sjbar.a0*Fj.a1*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum)))
  term5.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                            -t(apply(RMTLj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  term6.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                             -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  ueif.a1 <- 1/mean(bw*(a==1))*(term2.a1+term3.a1+term4.a1+term5.a1+term6.a1)+1/mean(tilt)*term1.a1

  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

#' @title Separable Indirect Effect for RMTLJ with A*=0
#'
#' @description
#' Computes efficient influence function for separable indirect effect of RMTL
#' with reference treatment A*=0 in competing risks mediation analysis.
#'
#' @param id (`numeric`) \\cr
#' Individual identifiers for observations.
#'
#' @param a (`numeric`) \\cr
#' Binary treatment indicator (0 or 1).
#'
#' @param time (`numeric`) \\cr
#' Observed time (minimum of event time and censoring time).
#'
#' @param event (`numeric`) \\cr
#' Event indicator (0 = censored, j = cause j event).
#'
#' @param bw (`numeric`) \\cr
#' Balancing weights (inverse propensity weights).
#'
#' @param tilt (`numeric`) \\cr
#' Tilting function for stabilization.
#'
#' @param G.a0, G.a1 (`matrix`) \\cr
#' Censoring survival function G(t|X,A=0) and G(t|X,A=1).
#'
#' @param S.a0, S.a1 (`matrix`) \\cr
#' Overall event-free survival function S(t|X,A=0) and S(t|X,A=1).
#'
#' @param Sj.a0, Sj.a1 (`matrix`) \\cr
#' Cause-specific survival for cause j: Sj(t|X,A=0) and Sj(t|X,A=1).
#'
#' @param Sjbar.a0, Sjbar.a1 (`matrix`) \\cr
#' Complementary cause survival (all causes except j): Sjbar(t|X,A=0) and Sjbar(t|X,A=1).
#'
#' @param freq.time (`numeric(1)`) \\cr
#' Time grid frequency. If NULL, uses unique event times.
#'
#' @param admin.cens (`numeric(1)`) \\cr
#' Administrative censoring time (maximum follow-up).
#'
#' @return List with elements:
#' \describe{
#'   \item{ueif.a1}{Efficient influence function for A=1 under separable indirect effect}
#'   \item{ueif.a0}{Efficient influence function for A=0 under separable indirect effect}
#' }
#'
#' @details
#' This function implements the separable indirect effect decomposition where the
#' treatment affects complementary causes while holding the effect on cause j at
#' reference level A*=0. The counterfactual RMTL combines Sj under A=0 with Sjbar
#' under the treatment level.
#'
#' The efficient influence function has 6 terms, symmetric to the A*=1 version but
#' with different complementary cause weighting for reference level A*=0.
#'
#' @keywords internal
double_ml_rmtlj_sep_indirect_astar0 <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))

  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Sjbar.a0 <- t(na.locf(t(ifelse(Sjbar.a0 < 1e-3, 1e-3, Sjbar.a0))))
  Sjbar.a1 <- t(na.locf(t(ifelse(Sjbar.a1 < 1e-3, 1e-3, Sjbar.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fjbar.dHazard.a0 <- t(apply(cbind(0, -log(Sjbar.a0)), 1, diff))
  Fjbar.dHazard.a1 <- t(apply(cbind(0, -log(Sjbar.a1)), 1, diff))

  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  Fj.a1.a0 <- t(apply(cbind(1, winsorize_values(Sj.a1*Sjbar.a0, lower = 1e-3, upper = 1)[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum)) # not Fj.dHazard.a1!

  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1.a0 <- t(apply(Fj.a1.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))

  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))

  # a=0
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1.a0
  term2.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply((dNjt-Yt*Fj.dHazard.a0)/(G.a0)*(1-Sjbar.a0/Sjbar.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNjt-Yt*Fj.dHazard.a0)/(G.a0)*(1-Sjbar.a0/Sjbar.a0), 1, cumsum)))
  term3.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a0*t(apply((dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a0/Sjbar.a0), 1, cumsum))
                                                             -t(apply(RMTLj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a0/Sjbar.a0), 1, cumsum)))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a0/Sjbar.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(1-Sjbar.a0/Sjbar.a0), 1, cumsum)))
  term5.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                             -t(apply(RMTLj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  term6.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  ueif.a0 <- 1/mean(bw*(a==0))*(term2.a0+term3.a0+term4.a0+term5.a0+term6.a0)+1/mean(tilt)*term1.a0

  # a=1
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1
  term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply((dNjt-Yt*Fj.dHazard.a1)/(G.a1)*(1-Sjbar.a0/Sjbar.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNjt-Yt*Fj.dHazard.a1)/(G.a1)*(1-Sjbar.a0/Sjbar.a1), 1, cumsum)))
  term3.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a1*t(apply((dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a0/Sjbar.a1), 1, cumsum))
                                                             -t(apply(RMTLj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a0/Sjbar.a1), 1, cumsum)))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a0/Sjbar.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(1-Sjbar.a0/Sjbar.a1), 1, cumsum)))
  term5.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(RMTLj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                             -t(apply(RMTLj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  term6.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                            -t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  ueif.a1 <- 1/mean(bw*(a==1))*(term2.a1+term3.a1+term4.a1+term5.a1+term6.a1)+1/mean(tilt)*term1.a1

  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

#' @title Double Machine Learning for Restricted Mean Time Lost in Competing Risks
#'
#' @description
#' R6 class for double machine learning estimation of weighted average treatment
#' effects (WATE) on restricted mean time lost (RMTL) to a specific cause in
#' competing risks settings.
#'
#' @details
#' This class implements double machine learning estimation for the restricted
#' mean time lost (RMTL) to cause j before time tau in competing risks analysis.
#' RMTL measures the expected time lost to a specific cause before the restriction
#' time tau.
#'
#' The estimator uses efficient influence functions (UEIF) with cross-fitting to
#' obtain doubly robust estimates of the average treatment effect:
#' \deqn{ATE(\tau) = E[RMTL_j(\tau)|A=1] - E[RMTL_j(\tau)|A=0]}
#'
#' Nuisance functions estimated via machine learning:
#' * Overall event survival: S(t|X,A)
#' * Cause-specific survival: Sj(t|X,A)
#' * Complementary cause survival: Sjbar(t|X,A)
#' * Censoring survival: G(t|X,A)
#' * Treatment propensity: P(A=1|X)
#'
#' @examples
#' \dontrun{
#' library(DoubleML)
#' library(SuperLearner)
#' library(survSuperLearner)
#'
#' # Generate competing risks data
#' set.seed(123)
#' data <- make_competing_data(n_obs = 2000, admin_cens = 10)
#'
#' # Define learner libraries
#' ml_g_surv <- c("survSL.km", "survSL.coxph")
#' ml_g_cens <- c("survSL.km", "survSL.coxph")
#' ml_m <- c("SL.mean", "SL.glm")
#'
#' # Estimate RMTL treatment effect
#' dml_rmtlj <- DoubleMLRMTLJ$new(
#'   data,
#'   ml_g_surv = ml_g_surv,
#'   ml_g_surv_j = ml_g_surv,
#'   ml_g_surv_jbar = ml_g_surv,
#'   ml_g_cens = ml_g_cens,
#'   ml_m = ml_m,
#'   admin_cens = 10,
#'   time_col = "time",
#'   event_col = "event",
#'   tau = 4,
#'   n_folds = 2
#' )
#'
#' dml_rmtlj$fit()
#' print(dml_rmtlj$coef)
#' print(dml_rmtlj$se)
#' }
#'
#' @export
DoubleMLRMTLJ = R6Class("DoubleMLRMTLJ",
  inherit = DoubleMLCompeting,
  public = list(
    #' @description
    #' Creates a new instance of this R6 class.
    #'
    #' @param data (`DoubleMLData`) \cr
    #' The `DoubleMLData` object providing the competing risks data.
    #'
    #' @param ml_g_surv (character vector) \cr
    #' survSuperLearner library for overall event survival S(t|X,A).
    #'
    #' @param ml_g_surv_j (character vector) \cr
    #' survSuperLearner library for cause-specific survival Sj(t|X,A).
    #'
    #' @param ml_g_surv_jbar (character vector) \cr
    #' survSuperLearner library for complementary cause survival Sjbar(t|X,A).
    #'
    #' @param ml_g_cens (character vector) \cr
    #' survSuperLearner library for censoring survival G(t|X,A).
    #'
    #' @param ml_m (character vector) \cr
    #' SuperLearner library for treatment propensity P(A=1|X).
    #'
    #' @param admin_cens (`numeric(1)`) \cr
    #' Administrative censoring time.
    #'
    #' @param time_col (`character(1)`) \cr
    #' Name of the time column in the data.
    #'
    #' @param event_col (`character(1)`) \cr
    #' Name of the event column (0=censored, 1=cause 1, 2=cause 2, ...).
    #'
    #' @param cause (`integer(1)`) \cr
    #' Cause of interest (default is 1).
    #'
    #' @param tau (`numeric(1)` or `numeric(vector)`) \cr
    #' Evaluation time point(s) for RMTL estimation. Default is admin_cens.
    #'
    #' @param freq_time (`numeric(1)`, optional) \cr
    #' Time grid frequency. If NULL, uses unique event times.
    #'
    #' @param ml_framework (`character(1)`) \cr
    #' Machine learning framework. Only "SuperLearner" supported. Default is "SuperLearner".
    #'
    #' @param n_folds (`integer(1)`) \cr
    #' Number of folds for cross-fitting. Default is 2.
    #'
    #' @param n_rep (`integer(1)`) \cr
    #' Number of repetitions for sample splitting. Default is 1.
    #'
    #' @param score (`character(1)`) \cr
    #' Score function. Default is "ueif".
    #'
    #' @param dml_procedure (`character(1)`) \cr
    #' DML procedure. Default is "dml2".
    #'
    #' @param draw_sample_splitting (`logical(1)`) \cr
    #' Whether to draw sample splitting. Default is TRUE.
    #'
    #' @param apply_cross_fitting (`logical(1)`) \cr
    #' Whether to apply cross-fitting. Default is TRUE.
    initialize = function(data,
                          ml_g_surv,
                          ml_g_surv_j,
                          ml_g_surv_jbar,
                          ml_g_cens,
                          ml_m,
                          admin_cens,
                          time_col,
                          event_col,
                          cause = 1,
                          tau = NULL,
                          freq_time = NULL,
                          ml_framework = "SuperLearner",
                          n_folds = 2,
                          n_rep = 1,
                          score = "ueif",
                          dml_procedure = "dml2",
                          draw_sample_splitting = TRUE,
                          apply_cross_fitting = TRUE) {

      # Set tau to admin_cens if not specified
      private$tau_ <- if (is.null(tau)) admin_cens else tau

      # Validate tau values
      private$validate_tau(admin_cens)

      # Initialize parent DoubleMLCompeting class
      super$initialize(
        data = data,
        ml_g_surv = ml_g_surv,
        ml_g_surv_j = ml_g_surv_j,
        ml_g_surv_jbar = ml_g_surv_jbar,
        ml_g_cens = ml_g_cens,
        ml_m = ml_m,
        admin_cens = admin_cens,
        time_col = time_col,
        event_col = event_col,
        cause = cause,
        freq_time = freq_time,
        ml_framework = ml_framework,
        n_folds = n_folds,
        n_rep = n_rep,
        score = score,
        dml_procedure = dml_procedure,
        draw_sample_splitting = draw_sample_splitting,
        apply_cross_fitting = apply_cross_fitting
      )
    },

    #' @description
    #' Compute score elements for DML estimation.
    #'
    #' @param ... Additional parameters (not used).
    #' @return List with psi_a and psi_b for DML framework.
    score_elements = function(...) {
      # Check if nuisance estimation has been performed
      if (is.null(private$all_ueif_)) {
        stop("Nuisance estimation must be performed before computing scores. Call fit() first.")
      }

      # Extract UEIF matrices from cross-fitting
      ueif_a1 <- private$all_ueif_$ueif_a1
      ueif_a0 <- private$all_ueif_$ueif_a0
      time_grid <- private$all_ueif_$time_grid

      # For single tau value
      if (length(private$tau_) == 1) {
        # Extract tau-specific estimates using extraction function
        raw_results <- rmtlj_extract_average(
          train.ueif.a1 = ueif_a1,
          train.ueif.a0 = ueif_a0,
          estimation.ueif.a1 = ueif_a1,
          estimation.ueif.a0 = ueif_a0,
          train.time = time_grid,
          estimation.time = time_grid,
          tau = private$tau_
        )

        # Wrap in list with tau
        inference_results <- list(
          tau = private$tau_,
          estimate.a0 = raw_results$estimate.a0,
          estimate.a1 = raw_results$estimate.a1,
          estimate.diff = raw_results$estimate.diff,
          se.a0 = raw_results$estimate.a0.se,
          se.a1 = raw_results$estimate.a1.se,
          se.diff = raw_results$estimate.diff.se
        )
      } else {
        # Handle multiple tau values
        inference_list <- lapply(private$tau_, function(tau_i) {
          rmtlj_extract_average(
            train.ueif.a1 = ueif_a1,
            train.ueif.a0 = ueif_a0,
            estimation.ueif.a1 = ueif_a1,
            estimation.ueif.a0 = ueif_a0,
            train.time = time_grid,
            estimation.time = time_grid,
            tau = tau_i
          )
        })

        # Combine results
        inference_results <- list(
          tau = private$tau_,
          estimate.a0 = sapply(inference_list, function(x) x$estimate.a0),
          estimate.a1 = sapply(inference_list, function(x) x$estimate.a1),
          estimate.diff = sapply(inference_list, function(x) x$estimate.diff),
          se.a0 = sapply(inference_list, function(x) x$estimate.a0.se),
          se.a1 = sapply(inference_list, function(x) x$estimate.a1.se),
          se.diff = sapply(inference_list, function(x) x$estimate.diff.se)
        )
      }

      # Store detailed RMTL estimates
      private$rmtlj_estimates_ <- inference_results

      # Extract individual-level influence functions at tau
      # Find time index closest to first tau
      tau_first <- private$tau_[1]
      tau_idx <- which.min(abs(time_grid - tau_first))

      # Individual influence functions at tau
      ueif_a1_at_tau <- ueif_a1[, tau_idx]
      ueif_a0_at_tau <- ueif_a0[, tau_idx]

      # Score for treatment effect: ψ(W) = ψ_1(W) - ψ_0(W)
      psi_b <- ueif_a1_at_tau - ueif_a0_at_tau

      # Jacobian component (for variance): ∂θ/∂ψ = -1 for linear functional
      psi_a <- rep(-1, nrow(ueif_a1))

      # Store coef and se for first tau (base class expects scalar)
      if (length(private$tau_) == 1) {
        private$coef_ <- inference_results$estimate.diff
        private$se_ <- inference_results$se.diff
      } else {
        private$coef_ <- inference_results$estimate.diff[1]
        private$se_ <- inference_results$se.diff[1]
      }

      return(list(psi_a = psi_a, psi_b = psi_b))
    },

    #' @description
    #' Compute separable direct effect with reference treatment A*=1.
    #'
    #' @details
    #' Estimates the direct effect of treatment on RMTL, holding the mediator at
    #' the distribution it would have under treatment A*=1. Requires the model
    #' to have been fitted first via fit().
    #'
    #' @return List with elements:
    #' \describe{
    #'   \item{tau}{Evaluation time point(s)}
    #'   \item{estimate.a0}{RMTL estimate under A=0, mediator at A*=1}
    #'   \item{estimate.a1}{RMTL estimate under A=1, mediator at A*=1}
    #'   \item{estimate.diff}{Direct effect (A=1 - A=0)}
    #'   \item{se.a0}{Standard error for A=0 estimate}
    #'   \item{se.a1}{Standard error for A=1 estimate}
    #'   \item{se.diff}{Standard error for direct effect}
    #' }
    fit_separable_direct_astar1 = function() {
      if (is.null(private$all_nuisance_)) {
        stop("Model must be fitted first. Call fit() before computing separable effects.")
      }

      # Use stored nuisance estimates and data
      n_obs <- self$data$n_obs
      time_grid <- private$all_ueif_$time_grid

      # Extract all stored data
      id <- 1:n_obs
      a <- private$all_nuisance_$a
      time <- private$all_nuisance_$time
      event <- private$all_nuisance_$event
      propensity <- private$all_nuisance_$propensity

      # Compute weights
      bw <- ifelse(a == 1, 1 / propensity, 1 / (1 - propensity))
      tilt <- rep(1, n_obs)

      # Call separable direct A*=1 function with stored nuisances
      ueif_result <- double_ml_rmtlj_sep_direct_astar1(
        id = id,
        a = a,
        time = time,
        event = event,
        bw = bw,
        tilt = tilt,
        G.a0 = private$all_nuisance_$G_a0,
        G.a1 = private$all_nuisance_$G_a1,
        S.a0 = private$all_nuisance_$S_a0,
        S.a1 = private$all_nuisance_$S_a1,
        Sj.a0 = private$all_nuisance_$Sj_a0,
        Sj.a1 = private$all_nuisance_$Sj_a1,
        Sjbar.a0 = private$all_nuisance_$Sjbar_a0,
        Sjbar.a1 = private$all_nuisance_$Sjbar_a1,
        freq.time = NULL,
        admin.cens = self$admin_cens
      )

      # Store UEIFs
      ueif_a1_sep <- ueif_result$ueif.a1
      ueif_a0_sep <- ueif_result$ueif.a0

      # Extract estimates at tau using inference function
      if (length(private$tau_) == 1) {
        raw_results <- rmtlj_extract_average(
          train.ueif.a1 = ueif_a1_sep,
          train.ueif.a0 = ueif_a0_sep,
          estimation.ueif.a1 = ueif_a1_sep,
          estimation.ueif.a0 = ueif_a0_sep,
          train.time = time_grid,
          estimation.time = time_grid,
          tau = private$tau_
        )

        results <- list(
          tau = private$tau_,
          estimate.a0 = raw_results$estimate.a0,
          estimate.a1 = raw_results$estimate.a1,
          estimate.diff = raw_results$estimate.diff,
          se.a0 = raw_results$estimate.a0.se,
          se.a1 = raw_results$estimate.a1.se,
          se.diff = raw_results$estimate.diff.se
        )
      } else {
        inference_list <- lapply(private$tau_, function(tau_i) {
          rmtlj_extract_average(
            train.ueif.a1 = ueif_a1_sep,
            train.ueif.a0 = ueif_a0_sep,
            estimation.ueif.a1 = ueif_a1_sep,
            estimation.ueif.a0 = ueif_a0_sep,
            train.time = time_grid,
            estimation.time = time_grid,
            tau = tau_i
          )
        })

        results <- list(
          tau = private$tau_,
          estimate.a0 = sapply(inference_list, function(x) x$estimate.a0),
          estimate.a1 = sapply(inference_list, function(x) x$estimate.a1),
          estimate.diff = sapply(inference_list, function(x) x$estimate.diff),
          se.a0 = sapply(inference_list, function(x) x$estimate.a0.se),
          se.a1 = sapply(inference_list, function(x) x$estimate.a1.se),
          se.diff = sapply(inference_list, function(x) x$estimate.diff.se)
        )
      }

      return(results)
    },

    #' @description
    #' Compute separable indirect effect with reference treatment A*=1.
    #'
    #' @details
    #' Estimates the indirect effect of treatment on RMTL through mediators,
    #' with treatment held at A*=1. Requires the model to have been fitted first.
    #'
    #' @return List with same structure as fit_separable_direct_astar1()
    fit_separable_indirect_astar1 = function() {
      if (is.null(private$all_nuisance_)) {
        stop("Model must be fitted first. Call fit() before computing separable effects.")
      }

      n_obs <- self$data$n_obs
      time_grid <- private$all_ueif_$time_grid
      id <- 1:n_obs
      a <- private$all_nuisance_$a
      time <- private$all_nuisance_$time
      event <- private$all_nuisance_$event
      propensity <- private$all_nuisance_$propensity
      bw <- ifelse(a == 1, 1 / propensity, 1 / (1 - propensity))
      tilt <- rep(1, n_obs)

      ueif_result <- double_ml_rmtlj_sep_indirect_astar1(
        id = id, a = a, time = time, event = event,
        bw = bw, tilt = tilt,
        G.a0 = private$all_nuisance_$G_a0, G.a1 = private$all_nuisance_$G_a1,
        S.a0 = private$all_nuisance_$S_a0, S.a1 = private$all_nuisance_$S_a1,
        Sj.a0 = private$all_nuisance_$Sj_a0, Sj.a1 = private$all_nuisance_$Sj_a1,
        Sjbar.a0 = private$all_nuisance_$Sjbar_a0, Sjbar.a1 = private$all_nuisance_$Sjbar_a1,
        freq.time = NULL, admin.cens = self$admin_cens
      )

      ueif_a1_sep <- ueif_result$ueif.a1
      ueif_a0_sep <- ueif_result$ueif.a0

      if (length(private$tau_) == 1) {
        raw_results <- rmtlj_extract_average(
          train.ueif.a1 = ueif_a1_sep, train.ueif.a0 = ueif_a0_sep,
          estimation.ueif.a1 = ueif_a1_sep, estimation.ueif.a0 = ueif_a0_sep,
          train.time = time_grid, estimation.time = time_grid,
          tau = private$tau_
        )

        results <- list(
          tau = private$tau_,
          estimate.a0 = raw_results$estimate.a0,
          estimate.a1 = raw_results$estimate.a1,
          estimate.diff = raw_results$estimate.diff,
          se.a0 = raw_results$estimate.a0.se,
          se.a1 = raw_results$estimate.a1.se,
          se.diff = raw_results$estimate.diff.se
        )
      } else {
        inference_list <- lapply(private$tau_, function(tau_i) {
          rmtlj_extract_average(
            train.ueif.a1 = ueif_a1_sep, train.ueif.a0 = ueif_a0_sep,
            estimation.ueif.a1 = ueif_a1_sep, estimation.ueif.a0 = ueif_a0_sep,
            train.time = time_grid, estimation.time = time_grid,
            tau = tau_i
          )
        })

        results <- list(
          tau = private$tau_,
          estimate.a0 = sapply(inference_list, function(x) x$estimate.a0),
          estimate.a1 = sapply(inference_list, function(x) x$estimate.a1),
          estimate.diff = sapply(inference_list, function(x) x$estimate.diff),
          se.a0 = sapply(inference_list, function(x) x$estimate.a0.se),
          se.a1 = sapply(inference_list, function(x) x$estimate.a1.se),
          se.diff = sapply(inference_list, function(x) x$estimate.diff.se)
        )
      }

      return(results)
    },

    #' @description
    #' Compute separable direct effect with reference treatment A*=0.
    #'
    #' @details
    #' Estimates the direct effect of treatment on RMTL, holding the mediator at
    #' the distribution it would have under treatment A*=0.
    #'
    #' @return List with same structure as fit_separable_direct_astar1()
    fit_separable_direct_astar0 = function() {
      if (is.null(private$all_nuisance_)) {
        stop("Model must be fitted first. Call fit() before computing separable effects.")
      }

      n_obs <- self$data$n_obs
      time_grid <- private$all_ueif_$time_grid
      id <- 1:n_obs
      a <- private$all_nuisance_$a
      time <- private$all_nuisance_$time
      event <- private$all_nuisance_$event
      propensity <- private$all_nuisance_$propensity
      bw <- ifelse(a == 1, 1 / propensity, 1 / (1 - propensity))
      tilt <- rep(1, n_obs)

      ueif_result <- double_ml_rmtlj_sep_direct_astar0(
        id = id, a = a, time = time, event = event,
        bw = bw, tilt = tilt,
        G.a0 = private$all_nuisance_$G_a0, G.a1 = private$all_nuisance_$G_a1,
        S.a0 = private$all_nuisance_$S_a0, S.a1 = private$all_nuisance_$S_a1,
        Sj.a0 = private$all_nuisance_$Sj_a0, Sj.a1 = private$all_nuisance_$Sj_a1,
        Sjbar.a0 = private$all_nuisance_$Sjbar_a0, Sjbar.a1 = private$all_nuisance_$Sjbar_a1,
        freq.time = NULL, admin.cens = self$admin_cens
      )

      ueif_a1_sep <- ueif_result$ueif.a1
      ueif_a0_sep <- ueif_result$ueif.a0

      if (length(private$tau_) == 1) {
        raw_results <- rmtlj_extract_average(
          train.ueif.a1 = ueif_a1_sep, train.ueif.a0 = ueif_a0_sep,
          estimation.ueif.a1 = ueif_a1_sep, estimation.ueif.a0 = ueif_a0_sep,
          train.time = time_grid, estimation.time = time_grid,
          tau = private$tau_
        )

        results <- list(
          tau = private$tau_,
          estimate.a0 = raw_results$estimate.a0,
          estimate.a1 = raw_results$estimate.a1,
          estimate.diff = raw_results$estimate.diff,
          se.a0 = raw_results$estimate.a0.se,
          se.a1 = raw_results$estimate.a1.se,
          se.diff = raw_results$estimate.diff.se
        )
      } else {
        inference_list <- lapply(private$tau_, function(tau_i) {
          rmtlj_extract_average(
            train.ueif.a1 = ueif_a1_sep, train.ueif.a0 = ueif_a0_sep,
            estimation.ueif.a1 = ueif_a1_sep, estimation.ueif.a0 = ueif_a0_sep,
            train.time = time_grid, estimation.time = time_grid,
            tau = tau_i
          )
        })

        results <- list(
          tau = private$tau_,
          estimate.a0 = sapply(inference_list, function(x) x$estimate.a0),
          estimate.a1 = sapply(inference_list, function(x) x$estimate.a1),
          estimate.diff = sapply(inference_list, function(x) x$estimate.diff),
          se.a0 = sapply(inference_list, function(x) x$estimate.a0.se),
          se.a1 = sapply(inference_list, function(x) x$estimate.a1.se),
          se.diff = sapply(inference_list, function(x) x$estimate.diff.se)
        )
      }

      return(results)
    },

    #' @description
    #' Compute separable indirect effect with reference treatment A*=0.
    #'
    #' @details
    #' Estimates the indirect effect of treatment on RMTL through mediators,
    #' with treatment held at A*=0.
    #'
    #' @return List with same structure as fit_separable_direct_astar1()
    fit_separable_indirect_astar0 = function() {
      if (is.null(private$all_nuisance_)) {
        stop("Model must be fitted first. Call fit() before computing separable effects.")
      }

      n_obs <- self$data$n_obs
      time_grid <- private$all_ueif_$time_grid
      id <- 1:n_obs
      a <- private$all_nuisance_$a
      time <- private$all_nuisance_$time
      event <- private$all_nuisance_$event
      propensity <- private$all_nuisance_$propensity
      bw <- ifelse(a == 1, 1 / propensity, 1 / (1 - propensity))
      tilt <- rep(1, n_obs)

      ueif_result <- double_ml_rmtlj_sep_indirect_astar0(
        id = id, a = a, time = time, event = event,
        bw = bw, tilt = tilt,
        G.a0 = private$all_nuisance_$G_a0, G.a1 = private$all_nuisance_$G_a1,
        S.a0 = private$all_nuisance_$S_a0, S.a1 = private$all_nuisance_$S_a1,
        Sj.a0 = private$all_nuisance_$Sj_a0, Sj.a1 = private$all_nuisance_$Sj_a1,
        Sjbar.a0 = private$all_nuisance_$Sjbar_a0, Sjbar.a1 = private$all_nuisance_$Sjbar_a1,
        freq.time = NULL, admin.cens = self$admin_cens
      )

      ueif_a1_sep <- ueif_result$ueif.a1
      ueif_a0_sep <- ueif_result$ueif.a0

      if (length(private$tau_) == 1) {
        raw_results <- rmtlj_extract_average(
          train.ueif.a1 = ueif_a1_sep, train.ueif.a0 = ueif_a0_sep,
          estimation.ueif.a1 = ueif_a1_sep, estimation.ueif.a0 = ueif_a0_sep,
          train.time = time_grid, estimation.time = time_grid,
          tau = private$tau_
        )

        results <- list(
          tau = private$tau_,
          estimate.a0 = raw_results$estimate.a0,
          estimate.a1 = raw_results$estimate.a1,
          estimate.diff = raw_results$estimate.diff,
          se.a0 = raw_results$estimate.a0.se,
          se.a1 = raw_results$estimate.a1.se,
          se.diff = raw_results$estimate.diff.se
        )
      } else {
        inference_list <- lapply(private$tau_, function(tau_i) {
          rmtlj_extract_average(
            train.ueif.a1 = ueif_a1_sep, train.ueif.a0 = ueif_a0_sep,
            estimation.ueif.a1 = ueif_a1_sep, estimation.ueif.a0 = ueif_a0_sep,
            train.time = time_grid, estimation.time = time_grid,
            tau = tau_i
          )
        })

        results <- list(
          tau = private$tau_,
          estimate.a0 = sapply(inference_list, function(x) x$estimate.a0),
          estimate.a1 = sapply(inference_list, function(x) x$estimate.a1),
          estimate.diff = sapply(inference_list, function(x) x$estimate.diff),
          se.a0 = sapply(inference_list, function(x) x$estimate.a0.se),
          se.a1 = sapply(inference_list, function(x) x$estimate.a1.se),
          se.diff = sapply(inference_list, function(x) x$estimate.diff.se)
        )
      }

      return(results)
    }
  ),

  active = list(
    #' @field coef Treatment effect estimate(s)
    coef = function() private$coef_,

    #' @field se Standard error of treatment effect estimate(s)
    se = function() private$se_,

    #' @field tau Evaluation time point(s)
    tau = function(value) {
      if (missing(value)) {
        return(private$tau_)
      } else {
        private$tau_ <- value
        private$validate_tau(private$admin_cens_)
      }
    },

    #' @field rmtlj_estimates All RMTL estimates and standard errors
    rmtlj_estimates = function() private$rmtlj_estimates_
  ),

  private = list(
    tau_ = NULL,
    coef_ = NULL,
    se_ = NULL,
    rmtlj_estimates_ = NULL,
    all_ueif_ = NULL,
    all_nuisance_ = NULL,  # Store nuisance estimates for separable effects

    validate_tau = function(admin_cens) {
      if (any(private$tau_ <= 0)) {
        stop("tau must be positive")
      }
      if (any(private$tau_ > admin_cens)) {
        stop("tau must be <= admin_cens")
      }
      invisible(NULL)
    },

    # Main nuisance estimation method called by base DoubleML class
    nuisance_est = function(smpls, ...) {
      # Extract full data
      full_data <- self$data$data
      n_obs <- self$data$n_obs

      # Get time grid
      time_grid <- self$get_time_grid()
      n_times <- length(time_grid)

      # Initialize prediction matrices for all observations
      S_a0_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      S_a1_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      Sj_a0_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      Sj_a1_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      Sjbar_a0_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      Sjbar_a1_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      G_a0_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      G_a1_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      propensity_pred <- rep(NA, n_obs)

      # Loop over folds for cross-fitting
      for (i_fold in seq_along(smpls$train_ids)) {
        train_ids <- smpls$train_ids[[i_fold]]
        test_ids <- smpls$test_ids[[i_fold]]

        # Get train/test data
        train_data <- full_data[train_ids, ]
        test_data <- full_data[test_ids, ]

        # Fit nuisance functions on train, predict on test
        ml_results <- private$ml_nuisance_est_superlearner(train_data, test_data, time_grid)

        # Store predictions for test fold
        S_a0_pred[test_ids, ] <- ml_results$S_a0
        S_a1_pred[test_ids, ] <- ml_results$S_a1
        Sj_a0_pred[test_ids, ] <- ml_results$Sj_a0
        Sj_a1_pred[test_ids, ] <- ml_results$Sj_a1
        Sjbar_a0_pred[test_ids, ] <- ml_results$Sjbar_a0
        Sjbar_a1_pred[test_ids, ] <- ml_results$Sjbar_a1
        G_a0_pred[test_ids, ] <- ml_results$G_a0
        G_a1_pred[test_ids, ] <- ml_results$G_a1
        propensity_pred[test_ids] <- ml_results$propensity
      }

      # Calculate balancing weights (inverse propensity weights)
      bw <- ifelse(full_data[[self$data$d_cols]] == 1,
                   1 / propensity_pred,
                   1 / (1 - propensity_pred))

      # Calculate tilting weights (for ATE, uniform weights)
      tilt <- rep(1, n_obs)

      # Compute UEIF using main RMTL function ON FULL SAMPLE with cross-fitted predictions
      ueif_result <- double_ml_rmtlj(
        id = 1:n_obs,
        a = full_data[[self$data$d_cols]],
        time = full_data[[self$time_col]],
        event = full_data[[self$event_col]],
        bw = bw,
        tilt = tilt,
        G.a0 = G_a0_pred,
        G.a1 = G_a1_pred,
        S.a0 = S_a0_pred,
        S.a1 = S_a1_pred,
        Sj.a0 = Sj_a0_pred,
        Sj.a1 = Sj_a1_pred,
        freq.time = private$freq_time_,
        admin.cens = self$admin_cens
      )

      # Store UEIF and nuisances for inference
      private$all_ueif_ <- list(
        ueif_a1 = ueif_result$ueif.a1,
        ueif_a0 = ueif_result$ueif.a0,
        time_grid = time_grid
      )

      # Store nuisance estimates for separable effects
      private$all_nuisance_ <- list(
        S_a0 = S_a0_pred,
        S_a1 = S_a1_pred,
        Sj_a0 = Sj_a0_pred,
        Sj_a1 = Sj_a1_pred,
        Sjbar_a0 = Sjbar_a0_pred,
        Sjbar_a1 = Sjbar_a1_pred,
        G_a0 = G_a0_pred,
        G_a1 = G_a1_pred,
        propensity = propensity_pred,
        a = full_data[[self$data$d_cols]],
        time = full_data[[self$time_col]],
        event = full_data[[self$event_col]]
      )

      # Compute scores via score_elements
      scores <- self$score_elements()

      # Return scores in DML format
      psi_elements <- list(psi_a = scores$psi_a, psi_b = scores$psi_b)
      return(psi_elements)
    }
  )
)
