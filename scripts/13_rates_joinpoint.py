"""Step 13: age-standardised rates and descriptive joinpoint analysis (protocol v1.0, sections 3, 5.3).

Age-standardised mortality rates (ASR, per 100,000, 2013 European Standard Population) by year for
Set L (WHO, 2003-2022, grid 0-14, 15-19, ..., 85+, common to all countries), its regions and countries
(MAL, BU, COMP, COMP_NARROW), and for Set W (malignant outcome, Eurostat 2011-2023 and UK from WHO;
grid 0-24, 25-29, ..., 85+).
Joinpoint (descriptive only): continuous piecewise log-linear model of ln(ASR), at most two joinpoints at
whole years, each segment at least four years long (end points included), weighted least squares with
weights 1/var(ln ASR) (Poisson variance of the ASR), number of joinpoints chosen by the BIC; APC per
segment and AAPC over the whole period with 95% intervals conditional on the joinpoint locations.
Primary period 2003-2022; S4: 2003-2019 (2020-2022 excluded).
Outputs: output/tables/asr_series.csv, joinpoint_segments.csv, joinpoint_aapc.csv
"""
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from config import OUTPUT, PROCESSED, ROOT, REGION
from lib.data import load_panel, collapse, esp_weights

TAB = OUTPUT / 'tables'
GRID_L = [0, 15] + list(range(20, 90, 5))
GRID_W = [0] + list(range(25, 90, 5))
OUTS = ['MAL', 'BU', 'COMP', 'COMP_NARROW']


def asr_table(df, grid, units, outs, by_year=True):
    w = esp_weights(grid)
    w = w / w.sum()
    rows = []
    for unit, members in units.items():
        s = df[df.geo.isin(members)].groupby(['year', 'age'], as_index=False)[outs + ['pop']].sum()
        for y, t in s.groupby('year'):
            t = t.set_index('age').reindex(w.index)
            rec = {'unit': unit, 'year': y, 'pop': t['pop'].sum()}
            for o in outs:
                r = t[o] / t['pop']
                rec[f'{o}_deaths'] = t[o].sum()
                rec[f'{o}_asr'] = float((w * r).sum() * 1e5)
                rec[f'{o}_asr_var'] = float((w ** 2 * t[o] / t['pop'] ** 2).sum() * 1e10)
            rows.append(rec)
    return pd.DataFrame(rows)


def design(t, knots):
    X = [np.ones_like(t), t]
    for k in knots:
        X.append(np.clip(t - k, 0, None))
    return np.column_stack(X)


def fit_joinpoint(years, y, var, max_jp=2, min_seg=4):
    t = np.asarray(years, float)
    w = 1 / np.asarray(var, float)
    n = len(t)
    best = None
    for k in range(max_jp + 1):
        cands = [()]
        if k:
            inner = [v for v in t[1:-1]]
            cands = [c for c in combinations(inner, k)
                     if all(seg >= min_seg for seg in np.diff(np.r_[t[0], c, t[-1]]) + 1)]
        for c in cands:
            X = design(t, c)
            W = np.sqrt(w)
            beta, *_ = np.linalg.lstsq(X * W[:, None], y * W, rcond=None)
            res = (y - X @ beta) * W
            sse = float(res @ res)
            p = 2 * k + 2
            bic = n * np.log(sse / n) + p * np.log(n)
            if best is None or bic < best['bic']:      # ties keep the model with fewer joinpoints
                best = {'k': k, 'knots': c, 'beta': beta, 'X': X, 'W': W, 'sse': sse, 'bic': bic}
    X, W, beta = best['X'], best['W'], best['beta']
    df_res = n - (2 * best['k'] + 2)
    s2 = best['sse'] / df_res if df_res > 0 else np.nan
    cov = s2 * np.linalg.inv((X * W[:, None]).T @ (X * W[:, None]))
    tq = stats.t.ppf(0.975, df_res) if df_res > 0 else np.nan
    bounds = np.r_[t[0], best['knots'], t[-1]]
    segs = []
    L = np.zeros((len(bounds) - 1, len(beta)))
    for i in range(len(bounds) - 1):
        L[i, 1] = 1
        L[i, 2:2 + i] = 1
        sl = float(L[i] @ beta)
        se = float(np.sqrt(L[i] @ cov @ L[i]))
        segs.append({'start': int(bounds[i]), 'end': int(bounds[i + 1]), 'apc': 100 * (np.exp(sl) - 1),
                     'apc_lo': 100 * (np.exp(sl - tq * se) - 1), 'apc_hi': 100 * (np.exp(sl + tq * se) - 1)})
    lens = np.diff(bounds)
    a = (lens[:, None] * L).sum(0) / lens.sum()
    sl, se = float(a @ beta), float(np.sqrt(a @ cov @ a))
    aapc = {'aapc': 100 * (np.exp(sl) - 1), 'aapc_lo': 100 * (np.exp(sl - tq * se) - 1),
            'aapc_hi': 100 * (np.exp(sl + tq * se) - 1), 'n_joinpoints': best['k'],
            'joinpoints': ' '.join(str(int(x)) for x in best['knots'])}
    return segs, aapc


def main():
    p = load_panel()
    cs = pd.read_csv(ROOT / 'protocol' / 'country_sets.csv')
    L = cs[cs.set_L].geo.tolist()
    rep = set(cs[cs.country_level_reporting].geo)
    q = collapse(p[p.geo.isin(L) & p.year.between(2003, 2022)], GRID_L, value_cols=OUTS + ['pop'])
    units = {'Set L': L, 'Set L excl. UK': [g for g in L if g != 'UK']}
    for r in sorted(set(REGION[g] for g in L)):
        units[r] = [g for g in L if REGION[g] == r]
    units.update({g: [g] for g in L})
    sL = asr_table(q, GRID_L, units, OUTS).assign(series='Set L (WHO)')
    # Set W, malignant outcome
    e = pd.read_csv(PROCESSED / 'panel_eurostat.csv.gz')
    e = e[e.sex.isin(['M', 'F'])].assign(sex=lambda x: x.sex.map({'M': 1, 'F': 2}))
    e = e[['geo', 'year', 'sex', 'age', 'MAL', 'pop_es']].rename(columns={'pop_es': 'pop'})
    uk = p[p.geo == 'UK'][['geo', 'year', 'sex', 'age', 'MAL', 'pop']]
    W = cs[cs.set_W].geo.tolist()
    ew = collapse(pd.concat([e, uk]), GRID_W, value_cols=['MAL', 'pop'])
    ew = ew[ew.year.between(2011, 2023) & ew.geo.isin(W)]
    complete = ew.groupby('geo').year.nunique()
    unitsW = {'Set W': W, 'Set W excl. UK': [g for g in W if g != 'UK']}
    sW = asr_table(ew, GRID_W, unitsW, ['MAL']).assign(series='Set W (Eurostat, UK WHO)')
    ny = ew.groupby('year').geo.nunique().rename('n_countries').reset_index()
    sW = sW.merge(ny, on='year')
    series = pd.concat([sL, sW], ignore_index=True)
    series.to_csv(TAB / 'asr_series.csv', index=False)

    segs, aap = [], []
    for (ser, unit), s in series.groupby(['series', 'unit']):
        if unit in L and unit not in rep:
            continue
        for o in (OUTS if ser.startswith('Set L') else ['MAL']):
            for label, yr in [('2003-2022', (2003, 2022)), ('S4 2003-2019', (2003, 2019)),
                              ('2011-2023', (2011, 2023)), ('2011-2019', (2011, 2019))]:
                if ser.startswith('Set L') != (yr[0] == 2003):
                    continue
                t = s[s.year.between(*yr)].sort_values('year')
                if ser.startswith('Set W'):
                    t = t[t.n_countries == t.n_countries.max()]
                if len(t) < 8 or (t[f'{o}_deaths'] == 0).any():
                    continue
                yv = np.log(t[f'{o}_asr'].to_numpy())
                var = t[f'{o}_asr_var'].to_numpy() / t[f'{o}_asr'].to_numpy() ** 2
                sg, ap = fit_joinpoint(t.year.to_numpy(), yv, var)
                meta = {'series': ser, 'unit': unit, 'outcome': o, 'period': label, 'n_years': len(t)}
                segs += [{**meta, **x} for x in sg]
                aap.append({**meta, **ap})
    pd.DataFrame(segs).to_csv(TAB / 'joinpoint_segments.csv', index=False)
    a = pd.DataFrame(aap)
    a.to_csv(TAB / 'joinpoint_aapc.csv', index=False)
    print(a[a.unit.isin(['Set L', 'Set L excl. UK', 'Set W'])].round(2).to_string(index=False))


if __name__ == '__main__':
    main()
