# Step 12: hierarchical negative binomial trend model (protocol v1.0, section 5.3).
#
# Data   : Set L, WHO panel (step 9), ages 0-14, 15-19, ..., 80-84, 85+; fitted on 2003-2019.
# Model  : deaths ~ age group x sex + calendar-year slopes by sex and broad age group (pooled, fixed)
#          + step terms at coding changes (country-specific level shifts, fixed)
#          + random intercepts and slopes for calendar year by country x sex x broad age group
#          (0-39, 40-64, 65-74, 75-84, 85+), unstructured 2 x 2 covariance; offset log person-years;
#          negative binomial (NB2) with log link. glmmTMB 1.1.8.
#          Age x sex interaction is included within the fixed effects for age group and sex.
# Step terms: documented and data-identified changes (data/processed/coding_changes.csv) in Set L
#          countries with 2004 <= year <= 2019 (a change in 2003, the first year, is not estimable).
# Outputs (output/tables): trend_pooled_slopes.csv, trend_variance.csv, trend_country_slopes.csv,
#          trend_steps.csv, trend_obs_exp.csv, trend_obs_exp_draws_summary.csv, trend_model_info.csv
# Expected deaths 2020-2023: 2,000 simulations; fixed and random effects drawn jointly from a
#          multivariate normal distribution with the joint precision matrix of the Laplace approximation
#          (TMB::sdreport, getJointPrecision = TRUE), deaths from the NB.
# Validation mode (additional analysis): FIT_END=<year> (2010, 2013, 2016) fits the same models to 2003-<year>
# and compares predictions for the following three years with observed deaths in countries without a coding
# change in those years; outputs carry the suffix _val<year> and models are not cached.
suppressPackageStartupMessages({library(glmmTMB); library(MASS); library(TMB); library(Matrix)})
fe <- as.integer(Sys.getenv('FIT_END', '2019'))
suf <- if (fe == 2019) '' else sprintf('_val%d', fe)
set.seed(20261029)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub('--file=', '', args[grep('--file=', args)])))
root <- dirname(here)
tab <- file.path(root, 'output', 'tables')
dir.create(tab, showWarnings = FALSE, recursive = TRUE)

p <- read.csv(gzfile(file.path(root, 'data', 'processed', 'panel_who.csv.gz')))
cs <- read.csv(file.path(root, 'protocol', 'country_sets.csv'))
L <- cs$geo[cs$set_L == 'True' | cs$set_L == TRUE]
p <- p[p$geo %in% L & p$year >= 2003 & p$year <= 2023 & !is.na(p$pop), ]
p$BU <- p$D32 + p$D33 + p$D42 + p$D43
p$COMP <- p$MAL + p$BU
p$agegrp <- ifelse(p$age < 15, 0, pmin(p$age, 85))
d <- aggregate(cbind(MAL, BU, COMP, pop) ~ geo + year + sex + agegrp, data = p, FUN = sum)
d$broad <- cut(d$agegrp, c(-1, 39, 64, 74, 84, 200), labels = c('0-39', '40-64', '65-74', '75-84', '85+'))
d$sexf <- factor(d$sex, 1:2, c('male', 'female'))
d$agef <- factor(d$agegrp)
d$year_c <- d$year - 2011
d$stratum <- interaction(d$geo, d$sexf, d$broad, drop = TRUE, sep = ':')

ch <- read.csv(file.path(root, 'data', 'processed', 'coding_changes.csv'))
ch_after <- unique(ch$geo[ch$geo %in% L & ch$year > fe & ch$year <= fe + 3])
ch <- ch[ch$geo %in% L & ch$year >= 2004 & ch$year <= fe, ]
step_names <- character(0)
for (i in seq_len(nrow(ch))) {
  nm <- sprintf('step_%s_%d', ch$geo[i], ch$year[i])
  d[[nm]] <- as.numeric(d$geo == ch$geo[i] & d$year >= ch$year[i])
  step_names <- c(step_names, nm)
}
fit_years <- d$year <= fe
fixed_rhs <- paste(c('agef * sexf', 'sexf:broad:year_c', step_names), collapse = ' + ')

fits <- list(); info <- list()
mdir <- file.path(root, 'output', 'models'); dir.create(mdir, showWarnings = FALSE)
for (o in c('MAL', 'BU', 'COMP')) {
  f <- as.formula(sprintf('%s ~ %s + offset(log(pop)) + (1 + year_c | stratum)', o, fixed_rhs))
  rds <- file.path(mdir, sprintf('trend_%s.rds', o))
  if (fe == 2019 && Sys.getenv('REUSE_MODELS') == '1' && file.exists(rds)) {   # development shortcut; default refits
    m <- readRDS(rds)
  } else {
    m <- glmmTMB(f, data = d[fit_years, ], family = nbinom2,
                 control = glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000)))
    if (fe == 2019) saveRDS(m, rds)
  }
  pd_ok <- isTRUE(m$sdr$pdHess)
  fits[[o]] <- m
  info[[o]] <- data.frame(outcome = o, n_obs = nobs(m), loglik = as.numeric(logLik(m)), AIC = AIC(m),
                          theta = sigma(m), pdHess = pd_ok, conv = m$fit$convergence,
                          n_strata = nlevels(droplevels(d$stratum[fit_years])), n_steps = length(step_names))
}
write.csv(do.call(rbind, info), file.path(tab, paste0('trend_model_info', suf, '.csv')), row.names = FALSE)

pooled <- list(); varc <- list(); cslopes <- list(); steps <- list()
for (o in names(fits)) {
  m <- fits[[o]]
  b <- fixef(m)$cond; V <- vcov(m)$cond
  sl <- grep('year_c', names(b), value = TRUE)
  for (s in sl) {
    est <- b[[s]]; se <- sqrt(V[s, s])
    parts <- regmatches(s, regexec('sexf(\\w+):broad([^:]+):year_c', s))[[1]]
    pooled[[length(pooled) + 1]] <- data.frame(outcome = o, sex = parts[2], broad_age = parts[3], beta = est, se = se,
      apc = 100 * (exp(est) - 1), apc_lo = 100 * (exp(est - 1.96 * se) - 1), apc_hi = 100 * (exp(est + 1.96 * se) - 1))
  }
  vc <- VarCorr(m)$cond$stratum
  varc[[o]] <- data.frame(outcome = o, sd_intercept = attr(vc, 'stddev')[1], sd_slope = attr(vc, 'stddev')[2],
                          corr = attr(vc, 'correlation')[1, 2])
  re <- ranef(m, condVar = TRUE)$cond$stratum
  cv <- attr(re, 'condVar')
  keys <- do.call(rbind, strsplit(rownames(re), ':'))
  for (k in seq_len(nrow(re))) {
    s <- sprintf('sexf%s:broad%s:year_c', keys[k, 2], keys[k, 3])
    est <- b[[s]] + re[k, 'year_c']
    se <- sqrt(V[s, s] + cv[2, 2, k])
    cslopes[[length(cslopes) + 1]] <- data.frame(outcome = o, geo = keys[k, 1], sex = keys[k, 2], broad_age = keys[k, 3],
      beta = est, se = se, apc = 100 * (exp(est) - 1), re_slope = re[k, 'year_c'], re_intercept = re[k, '(Intercept)'])
  }
  for (s in step_names) {
    est <- b[[s]]; se <- sqrt(V[s, s])
    steps[[length(steps) + 1]] <- data.frame(outcome = o, term = s, rr = exp(est), lo = exp(est - 1.96 * se),
                                             hi = exp(est + 1.96 * se), beta = est, se = se)
  }
}
write.csv(do.call(rbind, pooled), file.path(tab, paste0('trend_pooled_slopes', suf, '.csv')), row.names = FALSE)
write.csv(do.call(rbind, varc), file.path(tab, paste0('trend_variance', suf, '.csv')), row.names = FALSE)
write.csv(do.call(rbind, cslopes), file.path(tab, paste0('trend_country_slopes', suf, '.csv')), row.names = FALSE)
st <- do.call(rbind, steps)
st$geo <- sub('step_(\\w+)_(\\d+)', '\\1', st$term); st$year <- as.integer(sub('step_(\\w+)_(\\d+)', '\\2', st$term))
write.csv(st, file.path(tab, paste0('trend_steps', suf, '.csv')), row.names = FALSE)

# ---- observed versus expected, 2020-2023 -------------------------------------------------------
nsim <- 2000
new <- d[d$year > fe & d$year <= (if (fe == 2019) 2023 else min(fe + 3, 2019)), ]
if (fe != 2019) new <- new[!(new$geo %in% ch_after), ]
rows <- list()
for (o in names(fits)) {
  m <- fits[[o]]
  b <- fixef(m)$cond; V <- vcov(m)$cond
  X <- model.matrix(as.formula(paste('~', fixed_rhs)), data = new)
  X <- X[, names(b), drop = FALSE]
  re <- ranef(m, condVar = TRUE)$cond$stratum
  idx <- match(as.character(new$stratum), rownames(re))
  stopifnot(!anyNA(idx))
  # joint draws of fixed and random effects
  sdr <- sdreport(m$obj, getJointPrecision = TRUE)
  Q <- sdr$jointPrecision
  nmq <- rownames(Q)
  par <- m$obj$env$last.par.best
  stopifnot(all(names(par) == nmq))
  Lc <- Cholesky(Q, LDL = FALSE, perm = TRUE)
  z <- matrix(rnorm(nrow(Q) * nsim), nrow(Q))
  dr <- as.matrix(solve(Lc, solve(Lc, z, system = 'Lt'), system = 'Pt')) + par
  bi <- which(nmq == 'beta')
  stopifnot(length(bi) == length(b), max(abs(par[bi] - b)) < 1e-8)
  ui <- which(nmq == 'b')
  ns <- nrow(re)
  R0 <- dr[ui[seq(1, 2 * ns, 2)], ]; R1 <- dr[ui[seq(2, 2 * ns, 2)], ]
  stopifnot(max(abs(par[ui[seq(1, 2 * ns, 2)]] - re[, 1])) < 1e-8, max(abs(par[ui[seq(2, 2 * ns, 2)]] - re[, 2])) < 1e-8)
  eta <- X %*% dr[bi, ] + R0[idx, ] + R1[idx, ] * new$year_c + log(new$pop)
  mu <- exp(eta)
  th <- sigma(m)
  y <- matrix(rnbinom(length(mu), size = th, mu = mu), nrow = nrow(mu))
  mu_pt <- exp(as.numeric(X %*% b) + re[idx, 1] + re[idx, 2] * new$year_c + log(new$pop))
  grp <- interaction(new$geo, new$year, drop = TRUE, sep = '|')
  for (g in levels(grp)) {
    w <- which(grp == g)
    ys <- colSums(y[w, , drop = FALSE]); ms <- colSums(mu[w, , drop = FALSE])
    gy <- strsplit(g, '\\|')[[1]]
    rows[[length(rows) + 1]] <- data.frame(outcome = o, geo = gy[1], year = as.integer(gy[2]),
      observed = sum(new[[o]][w]), expected = sum(mu_pt[w]), exp_mean_lo = quantile(ms, .025), exp_mean_hi = quantile(ms, .975),
      pi_lo = quantile(ys, .025), pi_hi = quantile(ys, .975))
  }
  # aggregates over countries per year (each year: countries with observed data that year)
  for (yr in sort(unique(new$year))) {
    w <- which(new$year == yr)
    ys <- colSums(y[w, , drop = FALSE]); ms <- colSums(mu[w, , drop = FALSE])
    rows[[length(rows) + 1]] <- data.frame(outcome = o, geo = 'Set L', year = yr, observed = sum(new[[o]][w]),
      expected = sum(mu_pt[w]), exp_mean_lo = quantile(ms, .025), exp_mean_hi = quantile(ms, .975),
      pi_lo = quantile(ys, .025), pi_hi = quantile(ys, .975))
  }
  # 2020-2022 combined, countries observed in all three years
  w <- which(new$year <= fe + 3)
  ys <- colSums(y[w, , drop = FALSE]); ms <- colSums(mu[w, , drop = FALSE])
  rows[[length(rows) + 1]] <- data.frame(outcome = o, geo = 'Set L', year = as.integer(sprintf('%d%d', fe + 1, fe + 3)), observed = sum(new[[o]][w]),
    expected = sum(mu_pt[w]), exp_mean_lo = quantile(ms, .025), exp_mean_hi = quantile(ms, .975),
    pi_lo = quantile(ys, .025), pi_hi = quantile(ys, .975))
}
oe <- do.call(rbind, rows)
oe$ratio <- oe$observed / oe$expected
write.csv(oe, file.path(tab, paste0('trend_obs_exp', suf, '.csv')), row.names = FALSE)
cnt <- aggregate(geo ~ year, data = unique(new[, c('geo', 'year')]), FUN = length)
names(cnt)[2] <- 'n_countries'
write.csv(cnt, file.path(tab, paste0('trend_obs_exp_countries_per_year', suf, '.csv')), row.names = FALSE)
if (fe != 2019) write.csv(data.frame(excluded_change_in_window = ch_after), file.path(tab, paste0('trend_excluded', suf, '.csv')), row.names = FALSE)
sink(file.path(tab, paste0('trend_sessioninfo', suf, '.txt'))); print(sessionInfo()); sink()
print(do.call(rbind, info)); print(oe[oe$geo == 'Set L', ])
