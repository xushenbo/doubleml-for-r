# survival prob censoring martingale
double_ml_survival_prob <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  
  cens.martingale.integral.a0 <- t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a0 # cbind(1, G.a0[, 1:(ns-1)])
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*S.a0*(cens.martingale.integral.a0-1)
  ueif.a0 <- 1/mean(tilt)*(term1.a0+term2.a0+term3.a0)
  
  cens.martingale.integral.a1 <- t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a1 # cbind(1, G.a1[, 1:(ns-1)])
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*S.a1*(cens.martingale.integral.a1-1)
  ueif.a1 <- 1/mean(tilt)*(term1.a1+term2.a1+term3.a1)
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}


# rmst censoring martingale
double_ml_rmst <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  RMST.a0 <- t(apply(S.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMST.a1 <- t(apply(S.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Delta.tau.over.G.a0 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a0+matrix(event/rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Delta.tau.over.G.a1 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a1+matrix(event/rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Delta.tau.over.G.a0*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a0)/G.a0, 1, cumsum))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMST.a0*(t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))-1)
  term5.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  ueif.a0 <- 1/mean(tilt)*(term1.a0+term2.a0+term3.a0+term4.a0-term5.a0)
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Delta.tau.over.G.a1*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a1)/G.a1, 1, cumsum))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMST.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term5.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  ueif.a1 <- 1/mean(tilt)*(term1.a1+term2.a1+term3.a1+term4.a1-term5.a1)
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}
