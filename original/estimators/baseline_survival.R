# survival prob, outcome regression
or_survival_prob <- function(id, a, time, event, tilt, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a0
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  
  estimate.a0 <- 1/mean(tilt)*cummin(colMeans(term2.a0))
  estimate.a1 <- 1/mean(tilt)*cummin(colMeans(term2.a1))
  return(data.frame(s, estimate.a0, estimate.a1))
}


# survival prob, ipcw
ipcw_survival_prob <- function(id, a, time, event, bw, G.a0, G.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a0 # cbind(1, G.a0[, 1:(ns-1)])
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a1 # cbind(1, G.a1[, 1:(ns-1)])
  estimate.a0 <- 1/mean(bw*(a==0))*cummin(colMeans(term1.a0))
  estimate.a1 <- 1/mean(bw*(a==1))*cummin(colMeans(term1.a1))
  return(data.frame(s, estimate.a0, estimate.a1))
}

# survival prob, doubly robust
dr_survival_prob <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
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
  estimate.a0 <- cummin(1/mean(bw*(a==0))*colMeans(term1.a0+term3.a0)+1/mean(tilt)*colMeans(term2.a0))
  # 1/mean(bw*(a==0))*cumsum(cummin(colMeans(term1.a0+term3.a0))*ds)+1/mean(tilt)*cumsum(cummin(colMeans(term2.a0))*ds)
  
  cens.martingale.integral.a1 <- t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a1 # cbind(1, G.a1[, 1:(ns-1)])
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*S.a1*(cens.martingale.integral.a1-1)
  estimate.a1 <- cummin(1/mean(bw*(a==1))*colMeans(term1.a1+term3.a1)+1/mean(tilt)*colMeans(term2.a1))
  # 1/mean(bw*(a==1))*cumsum(cummin(colMeans(term1.a1+term3.a1))*ds)+1/mean(tilt)*cumsum(cummin(colMeans(term2.a1))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}


# rmst, outcome regression by integrating survival prob
or_rmst_integral <- function(id, a, time, event, tilt, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a0
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  
  estimate.a0 <- 1/mean(tilt)*cumsum(cummin(colMeans(term2.a0))*ds)
  estimate.a1 <- 1/mean(tilt)*cumsum(cummin(colMeans(term2.a1))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}

# rmst, ipcw by integrating survival prob
ipcw_rmst_integral <- function(id, a, time, event, bw, G.a0, G.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a0 # cbind(1, G.a0[, 1:(ns-1)])
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a1 # cbind(1, G.a1[, 1:(ns-1)])
  estimate.a0 <- 1/mean(bw*(a==0))*cumsum(cummin(colMeans(term1.a0))*ds)
  estimate.a1 <- 1/mean(bw*(a==1))*cumsum(cummin(colMeans(term1.a1))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}

# rmst, doubly robust by integrating survival prob
dr_rmst_integral <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
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
  estimate.a0 <- cumsum((cummin(1/mean(bw*(a==0))*colMeans(term1.a0+term3.a0)+1/mean(tilt)*colMeans(term2.a0)))*ds)
  # 1/mean(bw*(a==0))*cumsum(cummin(colMeans(term1.a0+term3.a0))*ds)+1/mean(tilt)*cumsum(cummin(colMeans(term2.a0))*ds)
  
  cens.martingale.integral.a1 <- t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a1 # cbind(1, G.a1[, 1:(ns-1)])
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*S.a1*(cens.martingale.integral.a1-1)
  estimate.a1 <- cumsum((cummin(1/mean(bw*(a==1))*colMeans(term1.a1+term3.a1)+1/mean(tilt)*colMeans(term2.a1)))*ds)
  # 1/mean(bw*(a==1))*cumsum(cummin(colMeans(term1.a1+term3.a1))*ds)+1/mean(tilt)*cumsum(cummin(colMeans(term2.a1))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}
