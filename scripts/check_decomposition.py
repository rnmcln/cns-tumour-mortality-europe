"""Independent check of the primary decomposition (step 11).

Recomputes the three-factor decomposition of the change in mean annual deaths between 2003-2005 and
2020-2022 (Das Gupta 1993, symmetric form) directly from data/processed/panel_who.csv.gz, without the
pipeline's decomposition library (lib/decomp.py) or data helpers (lib/data.py), and compares the
country and Set L totals with output/tables/decomp_primary.csv.

Deaths D = N * sum_i(s_i * r_i), with N the total person-years, s_i the share of person-years and r_i
the death rate in age-sex cell i. The symmetric three-factor components are
  population size  = (N1 - N0) * [ (S00 + S11) / 3 + (S01 + S10) / 6 ]
  age-sex structure = sum_i (s1_i - s0_i) * [ (N0 r0_i + N1 r1_i) / 3 + (N0 r1_i + N1 r0_i) / 6 ]
  rates            = sum_i (r1_i - r0_i) * [ (N0 s0_i + N1 s1_i) / 3 + (N0 s1_i + N1 s0_i) / 6 ]
where Sjk = sum_i sj_i * rk_i. The components add up to D1 - D0.

Usage: python3 check_decomposition.py   (exit status 1 if any difference exceeds 0.01 deaths)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GROUPS = {'MAL': ['MAL'], 'BU': ['D32', 'D33', 'D42', 'D43'], 'COMP': ['MAL', 'D32', 'D33', 'D42', 'D43']}
P0, P1 = range(2003, 2006), range(2020, 2023)
TOL = 0.01


def period_means(p, geo, years, cols):
    s = p[(p.geo == geo) & p.year.isin(years)]
    assert s.year.nunique() == len(years), (geo, list(years))
    g = s.groupby(['sex', 'age'])[cols + ['pop']].sum() / len(years)
    return g


def three_factor(d0, n0, d1, n1):
    N0, N1 = n0.sum(), n1.sum()
    s0, s1 = n0 / N0, n1 / N1
    r0, r1 = d0 / n0, d1 / n1
    S = lambda s, r: float((s * r).sum())  # noqa: E731
    size = (N1 - N0) * ((S(s0, r0) + S(s1, r1)) / 3 + (S(s0, r1) + S(s1, r0)) / 6)
    struct = float(((s1 - s0) * ((N0 * r0 + N1 * r1) / 3 + (N0 * r1 + N1 * r0) / 6)).sum())
    rates = float(((r1 - r0) * ((N0 * s0 + N1 * s1) / 3 + (N0 * s1 + N1 * s0) / 6)).sum())
    assert abs(size + struct + rates - (d1.sum() - d0.sum())) < 1e-6
    return size, struct, rates


def main():
    p = pd.read_csv(ROOT / 'data' / 'processed' / 'panel_who.csv.gz')
    sets = pd.read_csv(ROOT / 'protocol' / 'country_sets.csv')
    geos = sets[sets.set_L].geo.tolist()
    ref = pd.read_csv(ROOT / 'output' / 'tables' / 'decomp_primary.csv').set_index(['geo', 'outcome'])
    cols = sorted({c for g in GROUPS.values() for c in g})
    rows = []
    for geo in geos:
        a, b = period_means(p, geo, P0, cols), period_means(p, geo, P1, cols)
        assert a.index.equals(b.index), geo
        n0, n1 = a['pop'].to_numpy(float), b['pop'].to_numpy(float)
        for o, g in GROUPS.items():
            d0, d1 = a[g].sum(axis=1).to_numpy(float), b[g].sum(axis=1).to_numpy(float)
            size, struct, rates = three_factor(d0, n0, d1, n1)
            rows.append({'geo': geo, 'outcome': o, 'deaths_0': d0.sum(), 'deaths_1': d1.sum(),
                         'population_size': size, 'age_sex_structure': struct, 'rates': rates})
    res = pd.DataFrame(rows)
    tot = res.groupby('outcome')[['deaths_0', 'deaths_1', 'population_size', 'age_sex_structure', 'rates']].sum()
    tot = tot.assign(geo='Europe').reset_index()
    res = pd.concat([res, tot], ignore_index=True)
    res['behaviour_code'] = np.nan
    for geo in res.geo.unique():
        m = res.geo == geo
        res.loc[m & (res.outcome == 'MAL'), 'behaviour_code'] = (
            res.loc[m & (res.outcome == 'MAL'), 'rates'].iloc[0] - res.loc[m & (res.outcome == 'COMP'), 'rates'].iloc[0])
    keys = ['deaths_0', 'deaths_1', 'population_size', 'age_sex_structure', 'rates']
    diff = max(abs(getattr(r, k) - ref.loc[(r.geo, r.outcome), k]) for r in res.itertuples() for k in keys)
    coding = res[(res.geo == 'Europe') & (res.outcome == 'MAL')].behaviour_code.iloc[0]
    diff = max(diff, abs(coding - ref.loc[('Europe', 'CODING'), 'rates']))
    e = res[res.geo == 'Europe'].set_index('outcome')
    print('Set L totals, independent implementation (deaths per year):')
    print(e[keys].round(1).to_string())
    print(f'behaviour-code component: {coding:.1f}')
    print(f'maximum absolute difference from output/tables/decomp_primary.csv: {diff:.2e} deaths')
    sys.exit(0 if diff < TOL else 1)


if __name__ == '__main__':
    main()
