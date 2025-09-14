# cifj, outcome regression
or_cifj <- function(id, a, time, event, tilt, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a0
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1
  
  estimate.a0 <- 1/mean(tilt)*cummax(colMeans(term2.a0))
  estimate.a1 <- 1/mean(tilt)*cummax(colMeans(term2.a1))
  return(data.frame(s, estimate.a0, estimate.a1))
}

# cifj, ipcw
ipcw_cifj <- function(id, a, time, event, bw, G.a0, G.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Njt <- do.call(cbind, lapply(1:ns, function(u)  (time <= s[u])*(event == cause)))
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  estimate.a0 <- 1/mean(bw*(a==0))*cummax(colMeans(term1.a0))
  estimate.a1 <- 1/mean(bw*(a==1))*cummax(colMeans(term1.a1))
  return(data.frame(s, estimate.a0, estimate.a1))
}


# cifj, doubly robust
dr_cifj <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
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
  estimate.a0 <- cummax(1/mean(bw*(a==0))*colMeans(term1.a0+term3.a0-term4.a0)+1/mean(tilt)*colMeans(term2.a0))
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Fj.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNct - Yt*G.dHazard.a1)*Fj.a1/(S.a1*G.a1), 1, cumsum))
  estimate.a1 <- cummax(1/mean(bw*(a==1))*colMeans(term1.a1+term3.a1-term4.a1)+1/mean(tilt)*colMeans(term2.a1))
  return(data.frame(s, estimate.a0, estimate.a1))
}

# rmtlj, outcome regression by integrating cifj
or_rmtlj_integral <- function(id, a, time, event, tilt, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a0
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1
  
  estimate.a0 <- 1/mean(tilt)*cumsum(cummax(colMeans(term2.a0))*ds)
  estimate.a1 <- 1/mean(tilt)*cumsum(cummax(colMeans(term2.a1))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}

# rmtlj, ipcw by integrating cifj
ipcw_rmtlj_integral <- function(id, a, time, event, bw, G.a0, G.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Njt <- do.call(cbind, lapply(1:ns, function(u)  (time <= s[u])*(event == cause)))
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  estimate.a0 <- 1/mean(bw*(a==0))*cumsum(cummax(colMeans(term1.a0))*ds)
  estimate.a1 <- 1/mean(bw*(a==1))*cumsum(cummax(colMeans(term1.a1))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}

# rmtlj, dr by integrating cifj
dr_rmtlj_integral <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
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
  estimate.a0 <- cumsum((cummax(1/mean(bw*(a==0))*colMeans(term1.a0+term3.a0-term4.a0)+1/mean(tilt)*colMeans(term2.a0)))*ds)
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Fj.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNct - Yt*G.dHazard.a1)*Fj.a1/(S.a1*G.a1), 1, cumsum))
  estimate.a1 <- cumsum((cummax(1/mean(bw*(a==1))*colMeans(term1.a1+term3.a1-term4.a1)+1/mean(tilt)*colMeans(term2.a1)))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}

