source('app/model.R')
assert <- function(ok,message) if(!isTRUE(ok)) stop(message)
p <- default_config(); r <- simulate_model(p); d <- r$data
assert(abs(r$initial_slope-.196)<1e-12,'Baseline initial growth')
assert(abs(threshold(p)-50)<1e-12,'Baseline threshold')
assert(max(abs(rowSums(d[c('S','I','R')])-100))<1e-8,'SIR conservation')
assert(all(as.matrix(d[c('S','I','R')])>=0),'Nonnegative states')
peak <- which.max(d$I)
assert(peak>1 && peak<nrow(d),'Baseline epidemic peak')
assert(d$incidence[peak+1]>0 && d$infectious_entries[peak+1]<d$infectious_exits[peak+1], 'Transmission continues after peak')
p$beta <- .0015; low <- simulate_model(p)
assert(all(diff(low$data$I)<0),'Low transmission decline')
p <- default_config(); p$gamma <- .4
assert(simulate_model(p)$initial_slope<0 && threshold(p)==100,'Faster recovery')
p <- default_config(); p$I <- 0
assert(all(simulate_model(p)$data$I==0),'No spontaneous infection')
p <- default_config(); p$beta <- 0
r <- simulate_model(p)
assert(max(abs(r$data$I-exp(-p$gamma*r$data$time)))<1e-8,'Analytical recovery-only solution')
for(model in c('SIR','SIS','SEIR','SIRD')) {
 p <- default_config(); p$model <- model; p$m <- .03; p$v <- .2; p$E <- 5
 r <- simulate_model(p); fine <- simulate_model(p,max_step=.025)
 assert(max(abs(as.matrix(r$data[r$states])-as.matrix(fine$data[r$states])))<1e-5,paste(model,'step refinement'))
 total <- rowSums(r$data[r$states])
 assert(max(abs(total-total[1]))<1e-7,paste(model,'mass balance, including D'))
 if(model=='SIRD') assert(all(diff(r$data$D)>=0),'Cumulative deaths')
}
p <- default_config(); p$model <- 'SEIR'; p$E <- 0
assert(simulate_model(p)$initial_slope<0,'SEIR initial slope depends on E')
p$E <- 10
assert(simulate_model(p)$initial_slope>0,'SEIR exposed progression')
p <- default_config(); p$v <- 1
assert(identical(simulate_model(p)$data,simulate_model(default_config())$data),'Newborn vaccination has no effect without births')
p <- default_config(); p$S <- -1
assert(inherits(try(simulate_model(p),silent=TRUE),'try-error'),'Reject negative inputs')
p$S <- NA_real_
assert(inherits(try(simulate_model(p),silent=TRUE),'try-error'),'Reject missing inputs')
# Independent high-accuracy SciPy DOP853 fixtures generated in reference-check.py.
ref <- read.csv('tests/reference.csv')
for(model in unique(ref$model)) {
 p <- default_config(); p$model <- model
 r <- simulate_model(p)
 rows <- ref[ref$model==model,]
 for(i in seq_len(nrow(rows))) {
  row <- match(rows$time[i],r$data$time)
  assert(abs(r$data[[rows$state[i]]][row]-rows$value[i])<1e-5,paste(model,'independent reference'))
 }
}
parse('app/app.R')
cat('All model tests passed.\n')
