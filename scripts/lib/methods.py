"""Core methods: age standardisation, YLL, Shapley/Das Gupta decomposition, projections."""
from itertools import permutations

import numpy as np
import pandas as pd


def asmr(deaths, pop, std):
    """Age-standardised rate per 100,000 (direct method) and its Poisson standard error."""
    deaths, pop, w = np.asarray(deaths, float), np.asarray(pop, float), np.asarray(std, float)
    w = w / w.sum()
    rate = (w * deaths / pop).sum() * 1e5
    se = np.sqrt((w ** 2 * deaths / pop ** 2).sum()) * 1e5
    return rate, se


def shapley_three(f0, f1, combine):
    """Exact Shapley decomposition of combine(*factors) between state 0 and state 1.

    f0, f1: tuples of three factors (scalars or arrays); combine maps three factors to a scalar.
    Returns the three contributions; they sum exactly to combine(*f1) - combine(*f0).
    Equivalent to the Das Gupta (1993) symmetric decomposition for three factors.
    """
    k = 3
    contrib = np.zeros(k)
    perms = list(permutations(range(k)))
    for order in perms:
        state = list(f0)
        prev = combine(*state)
        for j in order:
            state[j] = f1[j]
            cur = combine(*state)
            contrib[j] += cur - prev
            prev = cur
    return contrib / len(perms)


def decompose_deaths(pop0, rate0, pop1, rate1, weight=None):
    """Decompose the change in deaths (or YLL if weight = life years per death by age).

    pop*, rate*: arrays by age group (population and age-specific death rate per person).
    Factors: population size (scalar), age structure (shares), age-specific rates.
    """
    pop0, pop1 = np.asarray(pop0, float), np.asarray(pop1, float)
    w = np.ones_like(pop0) if weight is None else np.asarray(weight, float)

    def combine(size, share, rate):
        return float((size * share * rate * w).sum())

    f0 = (pop0.sum(), pop0 / pop0.sum(), np.asarray(rate0, float))
    f1 = (pop1.sum(), pop1 / pop1.sum(), np.asarray(rate1, float))
    c = shapley_three(f0, f1, combine)
    return {'total_0': combine(*f0), 'total_1': combine(*f1),
            'population_size': c[0], 'age_structure': c[1], 'age_specific_rates': c[2]}


def attenuated_trend_projection(years, deaths, pop, proj_years, proj_pop, fit_from,
                                attenuation=(0.0, 0.25, 0.5, 0.75), period_len=5):
    """Age-specific log-linear Poisson trend, projected with Nordpred-style drift attenuation.

    For each age group the log rate is regressed on calendar year (fit_from onwards).
    The fitted slope is applied in full for the first projection period and reduced by the
    attenuation factors in successive periods of period_len years (cumulative on the log scale).
    deaths, pop: arrays [n_years, n_ages]; proj_pop: [n_proj_years, n_ages].
    Returns projected deaths [n_proj_years, n_ages] and the fitted slopes.
    """
    years = np.asarray(years)
    sel = years >= fit_from
    y, D, N = years[sel], np.asarray(deaths, float)[sel], np.asarray(pop, float)[sel]
    n_age = D.shape[1]
    slope = np.zeros(n_age)
    level = np.zeros(n_age)
    t = y - y.max()
    for a in range(n_age):
        b = _poisson_loglinear(t, D[:, a], N[:, a])
        level[a], slope[a] = b
    out = np.zeros((len(proj_years), n_age))
    for i, py in enumerate(proj_years):
        h = py - y.max()                         # years beyond the last observed year
        eff = 0.0
        for step in range(1, h + 1):              # accumulate attenuated drift year by year
            p = min((step - 1) // period_len, len(attenuation) - 1)
            eff += 1 - attenuation[p]
        out[i] = np.exp(level + slope * eff) * np.asarray(proj_pop, float)[i]
    return out, slope


def _poisson_loglinear(t, d, n, iters=50):
    """Poisson regression log(mu) = log(n) + b0 + b1 t by Newton-Raphson; returns (b0, b1)."""
    X = np.column_stack([np.ones_like(t, dtype=float), t.astype(float)])
    if d.sum() == 0:
        return np.log(1e-12), 0.0
    b = np.array([np.log(d.sum() / n.sum()), 0.0])
    for _ in range(iters):
        mu = n * np.exp(X @ b)
        g = X.T @ (d - mu)
        H = X.T @ (X * mu[:, None])
        step = np.linalg.solve(H, g)
        b = b + step
        if np.abs(step).max() < 1e-10:
            break
    return b[0], b[1]
