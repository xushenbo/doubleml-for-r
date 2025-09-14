## total effect
# cifj censoring martingale
double_ml_cifj <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
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
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  Njt <- do.call(cbind, lapply(1:ns, function(u)  (time <= s[u])*(event == cause)))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Fj.a0*(t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))-1)
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNct - Yt*G.dHazard.a0)*Fj.a0/(S.a0*G.a0), 1, cumsum))
  ueif.a0 <- 1/mean(tilt)*(term1.a0+term2.a0+term3.a0-term4.a0)
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Fj.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNct - Yt*G.dHazard.a1)*Fj.a1/(S.a1*G.a1), 1, cumsum))
  ueif.a1 <- 1/mean(tilt)*(term1.a1+term2.a1+term3.a1-term4.a1)
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

# rmtlj censoring martingale
double_ml_rmtlj <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
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
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  RMTLj.a0 <- t(apply(Fj.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMTLj.a1 <- t(apply(Fj.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  RMTLj.obs <- do.call(cbind, lapply(1:ns, function(u)  (s[u]-pmin(time, s[u]))*(event == cause)))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.obs/matrix(rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a0*(dNct - Yt*G.dHazard.a0)/(G.a0*S.a0), 1, cumsum))-
                                                              matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a0*(dNct - Yt*G.dHazard.a0)/(G.a0*S.a0), 1, cumsum)))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a0*(t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))-1)
  term5.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMTLj.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  ueif.a0 <- 1/mean(tilt)*(term1.a0+term2.a0+term3.a0+term4.a0-term5.a0)
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.obs/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*Fj.a1*(dNct - Yt*G.dHazard.a1)/(G.a1*S.a1), 1, cumsum))-
                                                              matrix(s,ncol=ns,nrow=n,byrow=TRUE)*t(apply(Fj.a1*(dNct - Yt*G.dHazard.a1)/(G.a1*S.a1), 1, cumsum)))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMTLj.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term5.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMTLj.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  ueif.a1 <- 1/mean(tilt)*(term1.a1+term2.a1+term3.a1+term4.a1-term5.a1)
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

## separable direct, astar=1
# cifj
double_ml_cifj_sep_direct_astar1 <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
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
  Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a0*Sjbar.a1))[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a1!
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))
  
  # a=0
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a0.a1
  term2.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(Sjbar.a1*(dNjt-Yt*Fj.dHazard.a0)/(Sjbar.a0*G.a0), 1, cumsum))
  term3.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a0*t(apply(Sjbar.a1*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum))
                                                             -t(apply(Fj.a0*Sjbar.a1*(dNt-Yt*S.dHazard.a0)/(Sjbar.a0*S.a0*G.a0), 1, cumsum)))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                            -t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  term5.a0 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                             -t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  ueif.a0 <- 1/mean(bw*(a==0))*(term2.a0+term3.a0+term4.a0)+1/mean(bw*(a==1))*term5.a0+1/mean(tilt)*term1.a0
  
  # a=1
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1
  term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNjt-Yt*Fj.dHazard.a1)/(G.a1), 1, cumsum))
  term3.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1*t(apply((dNt-Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                             -t(apply(Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  # term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(Sjbar.a1*(dNjt-Yt*Fj.dHazard.a1)/(Sjbar.a1*G.a1), 1, cumsum))
  # term3.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1*t(apply(Sjbar.a1*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum))
  #                                                            +t(apply(Fj.a1*Sjbar.a1*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum)))
  ueif.a1 <- 1/mean(bw*(a==1))*(term2.a1+term3.a1)+1/mean(tilt)*term1.a1
  
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

# rmtlj
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
  Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a0*Sjbar.a1))[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a1!
  
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

## separable indirect, astar=1
# cifj
double_ml_cifj_sep_indirect_astar1 <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
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
  Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a0*Sjbar.a1))[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a1!
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))
  
  # a=0
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a0
  term2.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNjt-Yt*Fj.dHazard.a0)/(G.a0), 1, cumsum))
  term3.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a0*t(apply((dNt-Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                             -t(apply(Fj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  ueif.a0 <- 1/mean(bw*(a==0))*(term2.a0+term3.a0)+1/mean(tilt)*term1.a0
  
  # a=1
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a0.a1
  term2.a1 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                            -t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  term3.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a0.a1*t(apply((dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                             -t(apply(Fj.a0.a1*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  term4.a1 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNjt-Yt*Fj.dHazard.a0)/(G.a0)*(Sjbar.a1/Sjbar.a0), 1, cumsum))
  term5.a1 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a0*t(apply((dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(Sjbar.a1/Sjbar.a0), 1, cumsum))
                                                             -t(apply(Fj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0)*(Sjbar.a1/Sjbar.a0), 1, cumsum)))
  ueif.a1 <- 1/mean(bw*(a==1))*(term3.a1)+1/mean(bw*(a==0))*(term2.a1+term4.a1+term5.a1)+1/mean(tilt)*term1.a1 # ueif.a1 <- 1/mean(bw*(a==1))*(term2.a1+term3.a1+term4.a1)+1/mean(tilt)*term1.a1
  
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

# rmtlj
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
  Fj.a0.a1 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a0*Sjbar.a1))[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum)) # not Fj.dHazard.a1!
  
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

## separable direct, astar=0
# cifj
double_ml_cifj_sep_direct_astar0 <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
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
  Fj.a1.a0 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a1*Sjbar.a0))[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum)) # not Fj.dHazard.a1!
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))
  
  # a=0
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a0
  term2.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNjt-Yt*Fj.dHazard.a0)/(G.a0), 1, cumsum))
  term3.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a0*t(apply((dNt-Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                             -t(apply(Fj.a0*(dNt-Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  ueif.a0 <- 1/mean(bw*(a==0))*(term2.a0+term3.a0)+1/mean(tilt)*term1.a0
  
  # a=1
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1.a0
  term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(Sjbar.a0*(dNjt-Yt*Fj.dHazard.a1)/(Sjbar.a1*G.a1), 1, cumsum))
  term3.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1*t(apply(Sjbar.a0*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum))
                                                             -t(apply(Fj.a1*Sjbar.a0*(dNt-Yt*S.dHazard.a1)/(Sjbar.a1*S.a1*G.a1), 1, cumsum)))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                            -t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  term5.a1 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                             -t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  ueif.a1 <- 1/mean(bw*(a==1))*(term2.a1+term3.a1+term4.a1)+1/mean(bw*(a==0))*term5.a1+1/mean(tilt)*term1.a1
  
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

# rmtlj
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
  Fj.a1.a0 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a1*Sjbar.a0))[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum)) # not Fj.dHazard.a1!
  
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

## separable indirect, astar=0
# cifj
double_ml_cifj_sep_indirect_astar0 <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sjbar.a0, Sj.a1, Sjbar.a1, freq.time=NULL, admin.cens)
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
  Fj.a1.a0 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a1*Sjbar.a0))[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum)) # not Fj.dHazard.a1!
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event != 0)))
  dNjt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  dNjbart <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event > 1)))
  
  # a=0
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1.a0
  term2.a0 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                            -t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  
  term3.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
                                                             -t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  term4.a0 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNjt-Yt*Fj.dHazard.a1)/(G.a1)*(Sjbar.a0/Sjbar.a1), 1, cumsum))
  term5.a0 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1*t(apply((dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(Sjbar.a0/Sjbar.a1), 1, cumsum))
                                                             -t(apply(Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(Sjbar.a0/Sjbar.a1), 1, cumsum)))
  ueif.a0 <- 1/mean(bw*(a==0))*(term3.a0)+1/mean(bw*(a==1))*(term2.a0+term4.a0+term5.a0)+1/mean(tilt)*term1.a0 # ueif.a1 <- 1/mean(bw*(a==1))*(term2.a1+term3.a1+term4.a1)+1/mean(tilt)*term1.a1
  
  # term2.a0 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNjt-Yt*Fj.dHazard.a1)/(G.a1)*(-Sjbar.a0/Sjbar.a1), 1, cumsum))
  # term3.a0 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1*t(apply((dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(-Sjbar.a0/Sjbar.a1), 1, cumsum))
  #                                                            -t(apply(Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1)*(-Sjbar.a0/Sjbar.a1), 1, cumsum)))
  # term4.a0 <- -matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  #                                                            -t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a0)/(S.a0*G.a0), 1, cumsum)))
  # term4.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1.a0*t(apply((dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  #                                                            -t(apply(Fj.a1.a0*(dNjbart-Yt*Fjbar.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  # ueif.a0 <- 1/mean(bw*(a==0))*(-term2.a0-term3.a0+term4.a0)-1/mean(bw*(a==1))*term4.a1+1/mean(tilt)*term1.a0
  
  # a=1
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1
  term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNjt-Yt*Fj.dHazard.a1)/(G.a1), 1, cumsum))
  term3.a1 <- -matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*(Fj.a1*t(apply((dNt-Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
                                                             -t(apply(Fj.a1*(dNt-Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum)))
  ueif.a1 <- 1/mean(bw*(a==1))*(term2.a1+term3.a1)+1/mean(tilt)*term1.a1
  
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

# rmtlj
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
  Fj.a1.a0 <- t(apply(cbind(1, Winsorize(minval = 1e-3, maxval = 1, (Sj.a1*Sjbar.a0))[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum)) # not Fj.dHazard.a1!
  
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


