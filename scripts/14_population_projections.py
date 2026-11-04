"""Step 14: projected person-years, 2020-2040 (protocol v1.0, section 5.4, populations).

- 30 study countries other than Serbia: EUROPOP2023 (proj_23np), 1 January populations by single year of
  age and sex; variants BSL (baseline), LMRT (lower mortality), HMIGR, LMIGR, NMIGR (higher, lower and no
  migration). Observed 1 January populations are used for 2020-2025 (demo_pjan). Because the observed
  2025 populations differ from EUROPOP2023 by up to 6% (audit/pop_projection_jump.csv), replacing observed
  years alone would create a discontinuity in 2025-2026; the projected age-sex-specific growth relative to
  2025 is therefore applied to the observed 2025 population (deviations table, entry 4, the method the protocol
  specifies for Serbia and the UK). Person-years = mean of 1 January of years t and t+1.
- Serbia and the United Kingdom: UN WPP 2024 medium variant (1 July, five-year groups to 100+);
  age-sex-specific growth factors WPP(t)/WPP(2024) applied to the last observed person-years (2024:
  Eurostat person-years for Serbia; WHO-reported ONS mid-year population for the UK).
- Observed person-years are used for 2020-2024 in all variants.
Ages are grouped to each country's grid (data/processed/audit/age_grids.csv) and, for the malignant
projection of Set W, to the Set W grid (0-24, 25-29, ..., top).
Outputs: data/processed/pop_projections.csv.gz (geo, variant, year, sex, age, py, grid),
         data/processed/audit/pop_projection_jump.csv (observed 2025 vs projected 2025 1 January totals).
"""
import numpy as np
import pandas as pd

from config import RAW, PROCESSED, COUNTRIES
from lib.data import load_panel

AUD = PROCESSED / 'audit'
VARIANTS = ['BSL', 'LMRT', 'HMIGR', 'LMIGR', 'NMIGR', 'EP25']   # EP25: EUROPOP2025 baseline (deviations table, entry 19)
YEARS = range(2020, 2041)


def grid_of(age, grid):
    g = np.asarray(grid)
    return g[np.searchsorted(g, np.asarray(age), side='right') - 1]


def single(a):
    if a == 'Y_LT1':
        return 0
    if a == 'Y_GE100':
        return 100
    if a.startswith('Y') and a[1:].isdigit():
        return int(a[1:])
    return None


def europop(grids):
    pr = pd.read_csv(RAW / 'eurostat' / 'proj_23np.csv.gz').rename(columns={'geo\\TIME_PERIOD': 'geo'})
    pr = pr[pr.sex.isin(['M', 'F']) & pr.projection.isin(VARIANTS)]
    p25 = pd.read_csv(RAW / 'eurostat' / 'proj_25np.csv.gz').rename(columns={'geo\\TIME_PERIOD': 'geo'})
    p25 = p25[p25.sex.isin(['M', 'F']) & (p25.projection == 'BSL')].assign(projection='EP25')
    pr = pd.concat([pr, p25[[c for c in pr.columns if c in p25.columns]]], ignore_index=True)
    pr['a1'] = pr.age.map(single)
    pr = pr[pr.a1.notna()]
    yrs = [str(y) for y in range(2022, 2042)]
    long = pr.melt(id_vars=['projection', 'geo', 'sex', 'a1'], value_vars=yrs, var_name='year', value_name='n')
    long['year'] = long.year.astype(int)
    ob = pd.read_csv(RAW / 'eurostat' / 'demo_pjan.csv.gz').rename(columns={'geo\\TIME_PERIOD': 'geo'})
    ob = ob[ob.sex.isin(['M', 'F']) & ob.geo.isin(COUNTRIES)]
    ob['a1'] = ob.age.map(lambda a: 0 if a == 'Y_LT1' else (int(a[1:]) if a.startswith('Y') and a[1:].isdigit() else (-1 if a == 'Y_OPEN' else None)))
    ob = ob[ob.a1.notna()].melt(id_vars=['geo', 'sex', 'a1'], value_vars=[str(y) for y in range(2020, 2026)],
                                 var_name='year', value_name='n_obs').dropna()
    ob['year'] = ob.year.astype(int)
    top = ob[ob.a1 >= 0].groupby(['geo', 'sex', 'year']).a1.max().rename('open') + 1
    ob = ob.merge(top.reset_index(), on=['geo', 'sex', 'year'])
    ob.loc[ob.a1 < 0, 'a1'] = ob.open
    jump, jump25 = [], []
    out = []
    for g in sorted(set(long.geo.unique()) | {'RS'}):
        if g not in grids:
            continue
        o = ob[ob.geo == g].copy()
        o['age'] = grid_of(o.a1, grids[g])
        jo = o.groupby(['sex', 'year', 'age']).n_obs.sum()
        if g == 'RS':                                  # observed person-years 2020-2024 only (WPP thereafter)
            jan = jo.unstack('year')
            py = (jan.loc[:, list(range(2020, 2025))].to_numpy() + jan.loc[:, list(range(2021, 2026))].to_numpy()) / 2
            py = pd.DataFrame(py, index=jan.index, columns=list(range(2020, 2025))).stack().rename('py').reset_index()
            py.columns = ['sex', 'age', 'year', 'py']
            out += [py.assign(geo=g, variant=v) for v in VARIANTS]
            continue
        for v in VARIANTS:
            s = long[(long.geo == g) & (long.projection == v)].copy()
            s['age'] = grid_of(s.a1, grids[g])
            jan = s.groupby(['sex', 'year', 'age']).n.sum()
            if v == 'EP25':
                jump25.append({'geo': g, 'projected_2025': jan.xs(2025, level='year').sum(),
                               'observed_2025': jo.xs(2025, level='year').sum()})
            if v == 'BSL':
                jump.append({'geo': g, 'projected_2025': jan.xs(2025, level='year').sum(),
                             'observed_2025': jo.xs(2025, level='year').sum()})
            # observed 1 January populations (2020-2025) are used; from 2026 the projected age-sex-specific
            # growth relative to 2025 is applied to the observed 2025 population (deviations table, entry 4)
            jp = jan.unstack('year')
            jobs = jo.unstack('year')
            fac = jp.loc[:, list(range(2026, 2042))].div(jp[2025], axis=0)
            jan = pd.concat([jobs.loc[:, list(range(2020, 2026))], fac.mul(jobs[2025], axis=0)], axis=1)
            assert jan.notna().all().all(), g
            py = (jan.loc[:, list(range(2020, 2041))].to_numpy() + jan.loc[:, list(range(2021, 2042))].to_numpy()) / 2
            py = pd.DataFrame(py, index=jan.index, columns=list(range(2020, 2041))).stack().rename('py').reset_index()
            py.columns = ['sex', 'age', 'year', 'py']
            out.append(py.assign(geo=g, variant=v))
    pd.DataFrame(jump).assign(rel_diff=lambda x: x.observed_2025 / x.projected_2025 - 1).to_csv(AUD / 'pop_projection_jump.csv', index=False)
    pd.DataFrame(jump25).assign(rel_diff=lambda x: x.observed_2025 / x.projected_2025 - 1).to_csv(AUD / 'pop_projection_jump_europop2025.csv', index=False)
    o = pd.concat(out, ignore_index=True)
    o['sex'] = o.sex.map({'M': 1, 'F': 2})
    return o


def wpp(grids, last_obs):
    w = pd.read_csv(RAW / 'unwpp' / 'wpp2024_srb_gbr_age5.csv.gz')
    w = w[w.Time.between(2024, 2040)]
    w['geo'] = w.ISO3_code.map({'SRB': 'RS', 'GBR': 'UK'})
    long = w.melt(id_vars=['geo', 'Time', 'AgeGrpStart'], value_vars=['PopMale', 'PopFemale'], var_name='s', value_name='n')
    long['sex'] = long.s.map({'PopMale': 1, 'PopFemale': 2})
    out = []
    for g in ['RS', 'UK']:
        s = long[long.geo == g].copy()
        s['age'] = grid_of(s.AgeGrpStart, grids[g])
        s = s.groupby(['sex', 'age', 'Time']).n.sum().unstack('Time')
        f = s.div(s[2024], axis=0)
        base = last_obs[last_obs.geo == g].set_index(['sex', 'age'])['pop']
        proj = f.mul(base, axis=0).loc[:, 2025:2040].stack().rename('py').reset_index().rename(columns={'Time': 'year'})
        for v in VARIANTS:
            out.append(proj.assign(geo=g, variant=v))
    return pd.concat(out, ignore_index=True)


def main():
    grids = {r.geo: [int(x) for x in r.grid.split()] for r in pd.read_csv(AUD / 'age_grids.csv').itertuples()}
    p = load_panel()
    ep = europop(grids)                                # 2020-2040 for 30 countries; 2020-2024 for Serbia
    uk = p[(p.geo == 'UK') & p.year.between(2020, 2024)][['geo', 'year', 'sex', 'age', 'pop']]
    last = pd.concat([uk[uk.year == 2024], ep[(ep.geo == 'RS') & (ep.year == 2024) & (ep.variant == 'BSL')]
                      .rename(columns={'py': 'pop'})[['geo', 'year', 'sex', 'age', 'pop']]])
    wp = wpp(grids, last)
    uk_v = pd.concat([uk.rename(columns={'pop': 'py'}).assign(variant=v) for v in VARIANTS])
    allp = pd.concat([ep, uk_v, wp], ignore_index=True)
    allp = allp[allp.year.isin(YEARS)].assign(grid='country')
    # check completeness: every country-variant-year has the same cells
    n = allp.groupby(['geo', 'variant', 'year']).size().groupby('geo').nunique()
    assert (n == 1).all(), n[n > 1]
    miss = allp.py.isna().sum()
    assert miss == 0, miss
    allp.to_csv(PROCESSED / 'pop_projections.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
    tot = allp[allp.variant == 'BSL'].groupby('year').py.sum()
    print('countries:', allp.geo.nunique(), '| total person-years (BSL):', tot.loc[[2020, 2024, 2025, 2030, 2040]].round(-3).to_dict())


if __name__ == '__main__':
    main()
