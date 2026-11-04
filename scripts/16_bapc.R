# Step 16: Bayesian age-period-cohort projections (protocol v1.0, section 5.4, model 5).
# BAPC 0.0.37 with R-INLA (Riebler and Held 2017): one model per country, outcome, origin and sex;
# ages 15 and over in five-year groups (country grid), annual periods (grid factor gf = 5),
# second-order random-walk priors for age, period and cohort (BAPC defaults) with an unstructured
# overdispersion term; predictions of age-specific counts (predictive mean and sd).
# Input : data/processed/bapc_jobs.csv.gz (step 15 prepare). Missing years in the history (Slovakia)
#         enter as missing counts. Ages 0-14 are projected by model 1 in step 15.
# Output: output/models/bapc_results.csv.gz (geo, outcome, origin, sex, age, year, mean, sd, pop)
suppressPackageStartupMessages({library(BAPC); library(INLA)})
INLA::inla.setOption(num.threads = '1:1')
set.seed(20261029)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub('--file=', '', args[grep('--file=', args)])))
root <- dirname(here)
j <- read.csv(gzfile(file.path(root, 'data', 'processed', 'bapc_jobs.csv.gz')))
keys <- unique(j[, c('geo', 'outcome', 'origin', 'sex')])
out <- vector('list', nrow(keys))
log <- data.frame()
for (k in seq_len(nrow(keys))) {
  s <- j[j$geo == keys$geo[k] & j$outcome == keys$outcome[k] & j$origin == keys$origin[k] & j$sex == keys$sex[k], ]
  ages <- sort(unique(s$age)); yrs <- sort(unique(s$year))
  Y <- matrix(NA_real_, length(yrs), length(ages), dimnames = list(yrs, ages))
  N <- Y
  Y[cbind(match(s$year, yrs), match(s$age, ages))] <- s$deaths
  N[cbind(match(s$year, yrs), match(s$age, ages))] <- s$pop
  npred <- sum(yrs > keys$origin[k])
  t0 <- Sys.time()
  res <- tryCatch({
    apc <- APCList(as.data.frame(Y), as.data.frame(N), gf = 5)
    r <- BAPC(apc, predict = list(npredict = npred, retro = TRUE), secondDiff = FALSE, stdweight = NULL,
              verbose = FALSE)
    pr <- agespec.proj(r)
    do.call(rbind, lapply(seq_along(ages), function(a) {
      x <- pr[[a]]
      data.frame(geo = keys$geo[k], outcome = keys$outcome[k], origin = keys$origin[k], sex = keys$sex[k],
                 age = ages[a], year = as.integer(rownames(x)), mean = x[, 'mean'], sd = x[, 'sd'],
                 pop = N[rownames(x), a])
    }))
  }, error = function(e) { message('FAILED ', paste(keys[k, ], collapse = ' '), ': ', conditionMessage(e)); NULL })
  if (!is.null(res)) res <- res[res$year > keys$origin[k], ]
  out[[k]] <- res
  log <- rbind(log, data.frame(keys[k, ], ok = !is.null(res), secs = as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
  if (k %% 25 == 0) cat(k, 'of', nrow(keys), '\n')
}
dir.create(file.path(root, 'output', 'models'), showWarnings = FALSE)
res <- do.call(rbind, out)
write.csv(res, gzfile(file.path(root, 'output', 'models', 'bapc_results.csv.gz')), row.names = FALSE)
write.csv(log, file.path(root, 'output', 'models', 'bapc_log.csv'), row.names = FALSE)
cat('fits:', nrow(log), 'failed:', sum(!log$ok), '\n')
sink(file.path(root, 'output', 'models', 'bapc_sessioninfo.txt')); print(sessionInfo()); sink()
