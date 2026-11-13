library(RTMB)
library(Matrix)

# sim setup params
r = 0.4
k = 1000
n_years = 30
q = 0.3
biom_pe = 0.1
srv_oe = 0.01
catch_oe = 0.2
u_ts = c(seq(0.1, 0.5, length.out = n_years / 2), rev(seq(0.1, 0.5, length.out = (n_years / 2) - 1)))
biom_devs = exp(rnorm(n_years-1, -biom_pe^2/2, biom_pe))

# mgmt quants
fmsy = r/2
msy = (r * k) / 4
bmsy = k / 2

# containers
biom = numeric(n_years)
catch = numeric(n_years-1)
obs_index = numeric(n_years)
obs_catch = numeric(n_years-1)

# initialize pop at carrying capacity
biom[1] = k

# loop through
for(t in 1:(n_years - 1)) {
  # popdy
  catch[t] = biom[t] * u_ts[t]
  biom[t+1] = (biom[t] + r*biom[t] * (1 - biom[t]/k) - catch[t]) * biom_devs[t]
  # generate observed
  obs_catch[t] = catch[t] * exp(rnorm(1,-catch_oe^2/2,catch_oe))
}
obs_index = q * biom * exp(rnorm(n_years,-srv_oe^2/2,srv_oe))


# make rtmb function
nll = function(pars) {

  RTMB::getAll(data, pars, warn = FALSE)

  # transform pars
  r = exp(ln_r)
  k = exp(ln_k)
  q = exp(ln_q)
  biom_pe = exp(ln_biom_pe)
  u_pe = exp(ln_u_pe)
  u_ts = plogis(logit_u_ts)

  # some derived quants
  msy = r * k / 4
  bmsy = k / 2
  umsy = r / 2

  # setup containers
  pred_biom = rep(0, n_years)
  pred_catch = rep(0, (n_years - 1))
  pred_index = rep(0, n_years)
  nll = 0

  # setup pop dy
  pred_biom[1] = k # first year

  # loop through after first year
  for(t in 1:(n_years - 1)) {
    pred_catch[t] = pred_biom[t] * u_ts[t]
    pred_biom[t+1] = (pred_biom[t] + r*pred_biom[t] * (1 - pred_biom[t]/k) - pred_catch[t]) * exp(ln_biom_devs[t])
  }
  pred_index = q * pred_biom

  # get nll
  nll = nll - sum(RTMB::dnorm(log(obs_catch), log(pred_catch) - catch_oe^2/2, catch_oe, TRUE))
  nll = nll - sum(RTMB::dnorm(log(obs_index), log(pred_index) - srv_oe^2/2, srv_oe, TRUE))
  nll = nll - sum(RTMB::dnorm(ln_biom_devs, -biom_pe^2/2, biom_pe, TRUE))
  nll = nll - sum(RTMB::dnorm(diff(logit_u_ts), 0, u_pe, TRUE))

  RTMB::REPORT(pred_biom)
  RTMB::REPORT(pred_catch)
  RTMB::REPORT(pred_index)
  RTMB::REPORT(u_ts)
  RTMB::ADREPORT(msy)
  RTMB::ADREPORT(bmsy)
  RTMB::ADREPORT(umsy)
  RTMB::ADREPORT(pred_biom)
  RTMB::ADREPORT(u_ts)

  return(nll)
}

data = list(
  obs_index = obs_index,
  obs_catch = obs_catch,
  catch_oe = catch_oe,
  srv_oe = srv_oe,
  n_years = n_years
)

pars = list(
  ln_r = log(r),
  ln_k = log(k),
  ln_q = log(q),
  ln_biom_pe = log(biom_pe),
  ln_u_pe = log(0.3),
  ln_biom_devs = log(biom_devs),
  logit_u_ts = qlogis(u_ts)
)

# map = list(ln_biom_pe = factor(NA))

obj = MakeADFun(nll, pars, random = c('ln_biom_devs', 'logit_u_ts'), map = NULL, silent = F)
optim = nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e4, eval.max = 1e4, rel.tol = 1e-10))
theta = optim$par
rep = obj$report(obj$env$last.par.best)
sdreport(obj)


