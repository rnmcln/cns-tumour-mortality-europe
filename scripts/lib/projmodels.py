"""Projection models (protocol v1.0, section 5.4). All models work on one country, sex-age cells as
columns: D_hist, N_hist [n_years, n_cells] for the fitting years, N_target [n_target, n_cells].
Each returns (mean [n_target, n_cells], draws [n_sim, n_target, n_cells]) where draws are expected
counts under parameter uncertainty (no Poisson noise; Poisson noise is added by the caller).

M1 constant: age-sex-specific rates of the last five observed years up to the origin
   (pooled deaths / pooled person-years); parameter draws rate ~ Gamma(D, 1 / N) (mean D / N; cells
   without deaths are projected at zero).
M2 damped log-linear (Moller et al. 2003 damping): per cell, Poisson regression of deaths on calendar
   year (offset log person-years) over the ten years up to the origin; the slope is applied in full in
   the first five projected years and reduced by 25%, 50% and 75% in successive five-year periods
   (cumulative on the log scale). Cells with fewer than 10 deaths in the fitting window use M1
   (implementation rule, not stated in the protocol). Parameter draws from the asymptotic normal
   distribution of (level, slope).
M3 age-drift: per sex, Poisson regression with age-group effects and one linear drift in calendar year
   over the ten years up to the origin, projected without damping; age groups without deaths in the
   window are projected at zero. Parameter draws from the asymptotic normal distribution.
"""
import numpy as np

FIT_YEARS = 10
CONST_YEARS = 5
MIN_DEATHS_M2 = 10
DAMP = (0.0, 0.25, 0.5, 0.75)


def damped_years(h, damp=DAMP, period=5):
    eff = 0.0
    for step in range(1, h + 1):
        eff += 1 - damp[min((step - 1) // period, len(damp) - 1)]
    return eff


def m1_constant(D, N, years, origin, Nt, rng, nsim):
    sel = (years <= origin) & (years > origin - CONST_YEARS)
    d, n = D[sel].sum(0), N[sel].sum(0)
    rate = np.where(n > 0, d / n, 0.0)
    mean = rate * Nt
    shape = np.where(d > 0, d, 1.0)
    draws = rng.gamma(shape=shape, scale=1.0 / n, size=(nsim, len(d))) * (d > 0)
    draws = draws[:, None, :] * Nt[None]
    return mean, draws


def _poisson_fit(X, d, n, iters=100):
    """Poisson regression with offset log n by Newton-Raphson; returns beta and covariance."""
    b = np.zeros(X.shape[1])
    b[0] = np.log(max(d.sum(), 0.5) / n.sum())
    for _ in range(iters):
        mu = n * np.exp(X @ b)
        g = X.T @ (d - mu)
        H = X.T @ (X * mu[:, None])
        step = np.linalg.solve(H, g)
        b += step
        if np.abs(step).max() < 1e-10:
            break
    mu = n * np.exp(X @ b)
    cov = np.linalg.inv(X.T @ (X * mu[:, None]))
    return b, cov


def m2_damped(D, N, years, origin, Nt, target_years, rng, nsim):
    sel = (years <= origin) & (years > origin - FIT_YEARS)
    y = years[sel] - origin
    Dw, Nw = D[sel], N[sel]
    m1_mean, m1_draws = m1_constant(D, N, years, origin, Nt, rng, nsim)
    mean = m1_mean.copy()
    draws = m1_draws.copy()
    eff = np.array([damped_years(int(t - origin)) for t in target_years])
    X = np.column_stack([np.ones(len(y)), y.astype(float)])
    slopes = np.full(D.shape[1], np.nan)
    for c in range(D.shape[1]):
        if Dw[:, c].sum() < MIN_DEATHS_M2:
            continue
        b, cov = _poisson_fit(X, Dw[:, c], Nw[:, c])
        slopes[c] = b[1]
        mean[:, c] = np.exp(b[0] + b[1] * eff) * Nt[:, c]
        bd = rng.multivariate_normal(b, cov, size=nsim)
        draws[:, :, c] = np.exp(bd[:, [0]] + bd[:, [1]] * eff[None, :]) * Nt[None, :, c]
    return mean, draws, slopes


def m3_agedrift(D, N, years, origin, Nt, target_years, sex_of_cell, rng, nsim):
    sel = (years <= origin) & (years > origin - FIT_YEARS)
    y = (years[sel] - origin).astype(float)
    Dw, Nw = D[sel], N[sel]
    mean = np.zeros_like(Nt, dtype=float)
    draws = np.zeros((nsim,) + Nt.shape)
    h = (np.asarray(target_years) - origin).astype(float)
    drift = {}
    for s in np.unique(sex_of_cell):
        cols = np.where((sex_of_cell == s) & (Dw.sum(0) > 0))[0]
        k = len(cols)
        rows_d, rows_n, Xr = [], [], []
        for j, c in enumerate(cols):
            for i in range(len(y)):
                x = np.zeros(k + 1)
                x[j] = 1
                x[k] = y[i]
                Xr.append(x)
                rows_d.append(Dw[i, c])
                rows_n.append(Nw[i, c])
        X, d, n = np.array(Xr), np.array(rows_d, float), np.array(rows_n, float)
        b = np.zeros(k + 1)
        b[:k] = np.log(Dw[:, cols].sum(0) / Nw[:, cols].sum(0))
        for _ in range(100):
            mu = n * np.exp(X @ b)
            step = np.linalg.solve(X.T @ (X * mu[:, None]), X.T @ (d - mu))
            b += step
            if np.abs(step).max() < 1e-10:
                break
        mu = n * np.exp(X @ b)
        cov = np.linalg.inv(X.T @ (X * mu[:, None]))
        drift[int(s)] = b[k]
        bd = rng.multivariate_normal(b, cov, size=nsim)
        for j, c in enumerate(cols):
            mean[:, c] = np.exp(b[j] + b[k] * h) * Nt[:, c]
            draws[:, :, c] = np.exp(bd[:, [j]] + bd[:, [k]] * h[None, :]) * Nt[None, :, c]
    return mean, draws, drift
