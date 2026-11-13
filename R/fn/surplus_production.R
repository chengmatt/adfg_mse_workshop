library(RTMB)
library(Matrix)

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
    pred_catch[t] = pred_biom[t] * u_ts[t] # get predicteed catch
    pred_biom[t+1] = (pred_biom[t] + r*pred_biom[t] * (1 - pred_biom[t]/k) - pred_catch[t]) * exp(ln_biom_devs[t]) # get pop dy
  }

  # get preedicted idx
  pred_index = q * pred_biom

  # get nll
  nll = nll - sum(RTMB::dnorm(log(obs_catch), log(pred_catch) - catch_oe^2/2, catch_oe, TRUE))
  nll = nll - sum(RTMB::dnorm(log(obs_index), log(pred_index[1:length(obs_index)]) - srv_oe^2/2, srv_oe, TRUE))
  nll = nll - sum(RTMB::dnorm(ln_biom_devs, -biom_pe^2/2, biom_pe, TRUE))
  nll = nll - sum(RTMB::dnorm(diff(logit_u_ts), 0, u_pe, TRUE))

  RTMB::REPORT(pred_biom)
  RTMB::REPORT(pred_catch)
  RTMB::REPORT(pred_index)
  RTMB::REPORT(u_ts)
  RTMB::REPORT(bmsy)
  RTMB::REPORT(umsy)
  RTMB::ADREPORT(msy)
  RTMB::ADREPORT(bmsy)
  RTMB::ADREPORT(umsy)
  RTMB::ADREPORT(pred_biom)
  RTMB::ADREPORT(u_ts)

  return(nll)
}

# sim setup params
r = 0.4
k = 1000
n_years = 30
n_cond_years = 15
q = 0.3
biom_pe = 0.1
srv_oe = 0.01
catch_oe = 0.2
u_cond_ts = c(seq(0.1, 0.5, length.out = n_cond_years / 2), rev(seq(0.1, 0.5, length.out = (n_cond_years / 2) - 1))) # harvest rate
biom_devs = exp(rnorm(n_years-1, -biom_pe^2/2, biom_pe)) # pe devs

# mgmt quants
fmsy = r/2
msy = (r * k) / 4
bmsy = k / 2

# containers
biom = numeric(n_years)
catch = numeric(n_years-1)
obs_index = numeric(n_years-1)
obs_catch = numeric(n_years-1)
u_ts = numeric(n_years - 1)
est_biom = rep(NA, n_years) # assessment estimate of current biomass in each mgmt year
tac = rep(NA, n_years - 1)

# initialize pop at carrying capacity
biom[1] = k

# loop through pop dy
for(t in 1:(n_years - 1)) {

  if(t <= n_cond_years) u_ts[t] = u_cond_ts[t] # conditioned based on hisotircal f
  else { # condition based on assessment estimates

    # setup data list
    data = list(
      obs_index = obs_index[1:(t-1)],
      obs_catch = obs_catch[1:(t-1)],
      catch_oe = catch_oe,
      srv_oe = srv_oe,
      n_years = t
    )

    # par list
    pars = list(
      ln_r = log(r),
      ln_k = log(k),
      ln_q = log(q),
      ln_biom_pe = log(biom_pe),
      ln_u_pe = log(0.3),
      ln_biom_devs = rep(0, t-1),
      logit_u_ts = qlogis(rep(0.2, t-1))
    )

    obj = MakeADFun(nll, pars, random = c('ln_biom_devs', 'logit_u_ts'), map = NULL, silent = T)
    optim = nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e4, eval.max = 1e4, rel.tol = 1e-10))
    theta = optim$par
    rep = obj$report(obj$env$last.par.best)

    # hcr: u = umsy if B > = bmsy, ramps linearly to 0 below bmsy
    est_biom[t] = rep$pred_biom[t] # save stuff
    tac[t] = rep$umsy * min(1, est_biom[t] / rep$bmsy) * est_biom[t]
    u_ts[t] = min(tac[t] / biom[t], 0.9)
  }

  # project pop forward w/ catch
  catch[t] = biom[t] * u_ts[t]
  biom[t+1] = (biom[t] + r*biom[t] * (1 - biom[t]/k) - catch[t]) * biom_devs[t]

  # generate obs for this year
  obs_index[t] = q * biom[t] * exp(rnorm(1,-srv_oe^2/2,srv_oe))
  obs_catch[t] = catch[t] * exp(rnorm(1,-catch_oe^2/2,catch_oe))

}

sdrep = sdreport(obj, getJointPrecision = TRUE)

plot(biom, type = 'l', ylim = c(0, k))
points(est_biom, col = 'red', pch = 16)
abline(h = bmsy, lty = 2)
# plot(rep$pred_biom, biom[1:length(rep$pred_biom)])
# abline(0,1)

