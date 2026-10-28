# Original implementation of the compartment equations documented by EcoEvoApps.
# Counts in a fixed study area; density-dependent incidence beta*S*I.
model_states <- function(model) switch(model, SIR=c('S','I','R'),
  SIS=c('S','I'), SEIR=c('S','E','I','R'), SIRD=c('S','I','R','D'),
  stop('Unknown model.'))
default_config <- function() list(model='SIR', S=99, E=0, I=1, R=0,
  beta=0.004, gamma=0.20, m=0, v=0, a=0.25, mu=0.05, days=60)
canonical_config <- function(p) {
  model_states(p$model)
  if(p$model!='SEIR') p$E <- 0
  if(p$model=='SIS') { p$R <- 0; p$v <- 0 }
  if(p$model!='SIRD') p$mu <- 0
  p
}
validate_config <- function(p) {
  p <- canonical_config(p)
  keys <- c('S','E','I','R','beta','gamma','m','v','a','mu','days')
  if(any(!vapply(p[keys], function(x) is.numeric(x) && length(x)==1 &&
                is.finite(x) && x>=0, logical(1)))) stop('Enter finite, nonnegative numbers in every field.')
  if(p$v>1) stop('Vaccinated fraction must be between 0 and 1.')
  if(p$days<1 || p$days>180) stop('Choose a duration from 1 to 180 days.')
  n <- sum(unlist(p[c('S','E','I','R')]))
  if(n<=0 || n>2000) stop('Initial living population must be greater than 0 and at most 2,000.')
  if(p$beta>0.02 || max(p$gamma,p$m,p$a,p$mu)>2) stop('Use beta <= 0.02 and other rates <= 2 per day.')
  p
}
rhs <- function(y,p) {
  S <- y[['S']]; I <- y[['I']]
  N <- sum(y[names(y)!='D']); infection <- p$beta*S*I
  z <- y*0
  z['S'] <- p$m*N*(1-p$v)-p$m*S-infection
  z['I'] <- infection-(p$gamma+p$m+p$mu)*I
  if(p$model=='SIS') {
    z['S'] <- z['S']+p$gamma*I
  } else {
    z['R'] <- p$gamma*I-p$m*y[['R']]+p$m*p$v*N
  }
  if(p$model=='SEIR') {
    z['E'] <- infection-(p$a+p$m)*y[['E']]
    z['I'] <- p$a*y[['E']]-(p$gamma+p$m)*I
  }
  if(p$model=='SIRD') z['D'] <- p$mu*I
  z
}
simulate_model <- function(config, output_step=0.25, max_step=0.1) {
  p <- validate_config(config)
  stopifnot(is.finite(output_step), output_step>0, is.finite(max_step), max_step>0)
  states <- model_states(p$model)
  init <- c(S=p$S,E=p$E,I=p$I,R=p$R,D=0)[states]
  times <- sort(unique(c(seq(0,p$days,by=output_step),p$days)))
  out <- matrix(NA_real_,nrow=length(times),ncol=length(states),dimnames=list(NULL,states))
  out[1,] <- y <- init
  # Bound the fastest flow to keep each explicit RK4 substep small.
  rate <- max(p$beta*sum(init)+p$m,p$a+p$m,p$gamma+p$m+p$mu,1)
  hmax <- min(max_step,0.2/rate)
  for(j in 2:length(times)) {
    n <- ceiling((times[j]-times[j-1])/hmax); h <- (times[j]-times[j-1])/n
    for(k in seq_len(n)) {
      k1 <- rhs(y,p); k2 <- rhs(y+h*k1/2,p)
      k3 <- rhs(y+h*k2/2,p); k4 <- rhs(y+h*k3,p)
      y <- y+h*(k1+2*k2+2*k3+k4)/6
      if(any(!is.finite(y)) || any(y < -1e-8)) stop('Numerical integration failed; reduce the rates.')
    }
    out[j,] <- y
  }
  dat <- data.frame(time=times,out,check.names=FALSE)
  dat$incidence <- p$beta*dat$S*dat$I
  dat$recoveries <- p$gamma*dat$I
  dat$infectious_exits <- (p$gamma+p$m+p$mu)*dat$I
  dat$infectious_entries <- if(p$model=='SEIR') p$a*dat$E else dat$incidence
  list(config=p,data=dat,states=states,initial_slope=unname(rhs(init,p)['I']))
}
threshold <- function(p) {
  p <- canonical_config(p)
  if(p$beta==0) return(Inf)
  if(p$model=='SEIR') {
    if(p$a==0) return(Inf)
    return((p$gamma+p$m)*(p$a+p$m)/(p$beta*p$a))
  }
  (p$gamma+p$m+p$mu)/p$beta
}
