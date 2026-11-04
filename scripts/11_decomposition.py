"""Step 11: decomposition of the change in recorded CNS tumour deaths (protocol v1.0, section 5.1).

Primary analysis: Set L, mean annual deaths 2003-2005 against 2020-2022, by country, outcome
(MAL = C70-C72; COMP = C70-C72, D32-D33, D42-D43; COMP_NARROW = C70-C72, D33, D43; BU), into population
size, age-sex structure and age-sex-specific rates (lib/decomp.py). Coding component = rate component of
MAL minus rate component of COMP; residual rate component = rate component of COMP. Country components
are summed to UN M49 regions and Europe (with and without the UK).
Uncertainty: Poisson bootstrap of the death counts of every country-year-sex-age cell and elementary
group (2,000 replicates, seed config.SEED), percentile 95% intervals; systematic: sensitivity analyses.
Also: chained annual decomposition (2003-2022, Set L), and Set W (MAL, Eurostat, 2011-2013 against
2021-2023; UK from WHO).

Sensitivity analyses (protocol 5.6) computed here: S1 narrow composite, S2 ages truncated at 85+,
S3 end period 2017-2019, S6 Set W with WHO counts where available, S7 WHO population denominators,
S8 exclusion of countries with a coding change between the periods, S9 complete-case subset,
S10 with and without the UK; additional: Germany 2021-2022 C70-C72 from Eurostat age-sex cells (deviations table, entry 2).
Outputs: output/tables/decomp_*.csv
"""
import numpy as np
import pandas as pd

from config import OUTPUT, PROCESSED, ROOT, REGION, SEED, N_BOOT, BASE_L, END_L, BASE_W, END_W
from lib.data import load_panel, collapse
from lib.decomp import decompose

TAB = OUTPUT / 'tables'
ELEM = ['MAL', 'D32', 'D33', 'D42', 'D43']
OUT_GROUPS = {'MAL': ['MAL'], 'COMP': ELEM, 'COMP_NARROW': ['MAL', 'D33', 'D43'], 'BU': ['D32', 'D33', 'D42', 'D43'],
              'D32': ['D32'], 'D33': ['D33'], 'D42': ['D42'], 'D43': ['D43']}
COMP_KEYS = ['deaths_0', 'deaths_1', 'population_size', 'age_sex_structure', 'rates']


def years(p):
    return list(range(p[0], p[1] + 1))


def cells(p, geo, yrs, popcol='pop'):
    """Arrays [n_years, n_cells] of elementary deaths and person-years for one country."""
    s = p[(p.geo == geo) & p.year.isin(yrs)].sort_values(['year', 'sex', 'age'])
    ny = len(yrs)
    assert s.year.nunique() == ny, (geo, yrs)
    nc = len(s) // ny
    d = {g: s[g].to_numpy(float).reshape(ny, nc) for g in ELEM}
    n = s[popcol].to_numpy(float).reshape(ny, nc)
    assert np.isfinite(n).all() and (n > 0).all(), (geo, yrs, popcol)
    key = s[['sex', 'age']].iloc[:nc].to_numpy()
    return d, n, key


def country_decomp(p, geo, base, end, popcol='pop', rng=None, nboot=0):
    d0, n0, k0 = cells(p, geo, years(base), popcol)
    d1, n1, k1 = cells(p, geo, years(end), popcol)
    assert (k0 == k1).all()
    N0, N1 = n0.mean(0), n1.mean(0)
    res = {}
    for o, grp in OUT_GROUPS.items():
        D0 = sum(d0[g] for g in grp).mean(0)
        D1 = sum(d1[g] for g in grp).mean(0)
        res[o] = {k: float(v) for k, v in decompose(D0, N0, D1, N1).items()}
    boot = None
    if nboot:
        sim0 = {g: rng.poisson(np.broadcast_to(d0[g], (nboot,) + d0[g].shape)).mean(1) for g in ELEM}
        sim1 = {g: rng.poisson(np.broadcast_to(d1[g], (nboot,) + d1[g].shape)).mean(1) for g in ELEM}
        boot = {o: decompose(sum(sim0[g] for g in grp), N0, sum(sim1[g] for g in grp), N1)
                for o, grp in OUT_GROUPS.items()}
    return res, boot


def rate_by_age(p, geos):
    """Cell-level components of the primary decomposition summed by sex and broad age group (Set L)."""
    rows = []
    for g in geos:
        d0, n0, key = cells(p, g, years(BASE_L))
        d1, n1, _ = cells(p, g, years(END_L))
        N0, N1 = n0.mean(0), n1.mean(0)
        for o, grp in OUT_GROUPS.items():
            r = decompose(sum(d0[x] for x in grp).mean(0), N0, sum(d1[x] for x in grp).mean(0), N1, by_cell=True)
            for i, (sx, age) in enumerate(key):
                rows.append({'geo': g, 'outcome': o, 'sex': int(sx), 'age': int(age),
                             **{k: float(v[i]) for k, v in r.items()}})
    c = pd.DataFrame(rows)
    c['broad'] = pd.cut(c.age, [-1, 39, 64, 74, 84, 200], labels=['0-39', '40-64', '65-74', '75-84', '85+']).astype(str)
    num = ['deaths_0', 'deaths_1', 'population_size', 'age_sex_structure', 'rates']
    t = c.groupby(['outcome', 'sex', 'broad'], as_index=False)[num].sum()
    cod = t[t.outcome == 'MAL'].set_index(['sex', 'broad']).rates - t[t.outcome == 'COMP'].set_index(['sex', 'broad']).rates
    t = pd.concat([t, cod.rename('rates').reset_index().assign(outcome='CODING')], ignore_index=True)
    return t


def tidy(res_by_geo, spec):
    rows = []
    for geo, res in res_by_geo.items():
        for o, r in res.items():
            rows.append({'spec': spec, 'geo': geo, 'outcome': o, **r})
        rows.append({'spec': spec, 'geo': geo, 'outcome': 'CODING', 'rates': res['MAL']['rates'] - res['COMP']['rates']})
    return pd.DataFrame(rows)


def aggregate(df, include_uk=True):
    d = df if include_uk else df[df.geo != 'UK']
    d = d.assign(region=d.geo.map(REGION))
    num = COMP_KEYS
    reg = d.groupby(['spec', 'region', 'outcome'], as_index=False)[num].sum(min_count=1).rename(columns={'region': 'geo'})
    eu = d.groupby(['spec', 'outcome'], as_index=False)[num].sum(min_count=1).assign(geo='Europe' if include_uk else 'Europe excl. UK')
    if not include_uk:
        reg['geo'] = reg.geo + ' excl. UK'
    return pd.concat([reg, eu], ignore_index=True)


def run_spec(p, geos, base, end, spec, popcol='pop', outcome_map=None):
    global OUT_GROUPS
    saved = OUT_GROUPS
    if outcome_map:
        OUT_GROUPS = outcome_map
    try:
        res = {g: country_decomp(p, g, base, end, popcol)[0] for g in geos}
    finally:
        OUT_GROUPS = saved
    df = tidy(res, spec)
    return pd.concat([df, aggregate(df.drop(columns=[]), True), aggregate(df, False)], ignore_index=True)


def primary_with_bootstrap(p, geos):
    rng = np.random.default_rng(SEED)
    point, boots = {}, {}
    for g in geos:                                   # fixed order -> reproducible streams
        point[g], boots[g] = country_decomp(p, g, BASE_L, END_L, rng=rng, nboot=N_BOOT)
    df = tidy(point, 'primary')
    allrows = [df, aggregate(df, True), aggregate(df, False)]
    point_tab = pd.concat(allrows, ignore_index=True)
    # bootstrap intervals for countries, regions, Europe (with/without UK)
    units = {g: [g] for g in geos}
    for reg in sorted(set(REGION[g] for g in geos)):
        units[reg] = [g for g in geos if REGION[g] == reg]
        units[reg + ' excl. UK'] = [g for g in geos if REGION[g] == reg and g != 'UK']
    units['Europe'] = list(geos)
    units['Europe excl. UK'] = [g for g in geos if g != 'UK']
    rows = []
    for u, members in units.items():
        if not members:
            continue
        for o in list(OUT_GROUPS) + ['CODING']:
            for k in ['population_size', 'age_sex_structure', 'rates', 'change']:
                if o == 'CODING':
                    if k != 'rates':
                        continue
                    v = sum(boots[g]['MAL']['rates'] - boots[g]['COMP']['rates'] for g in members)
                elif k == 'change':
                    v = sum(boots[g][o]['deaths_1'] - boots[g][o]['deaths_0'] for g in members)
                else:
                    v = sum(boots[g][o][k] for g in members)
                lo, hi = np.percentile(v, [2.5, 97.5])
                rows.append({'geo': u, 'outcome': o, 'component': k, 'boot_lo': lo, 'boot_hi': hi, 'boot_sd': v.std(ddof=1)})
    return point_tab, pd.DataFrame(rows)


def chained(p, geos, first=2003, last=2022):
    """Annual decompositions between consecutive observed years, accumulated from `first`.
    A step that spans a missing year (Slovakia 2011, 2015) is assigned to the later year."""
    rows = []
    for g in geos:
        obs = sorted(p[(p.geo == g) & p.year.between(first, last)].year.unique())
        cum = {o: np.zeros(3) for o in OUT_GROUPS}
        for o in OUT_GROUPS:
            rows.append({'geo': g, 'year': obs[0], 'outcome': o, 'population_size': 0.0, 'age_sex_structure': 0.0,
                         'rates': 0.0, 'step_from': obs[0]})
        for y0, y1 in zip(obs[:-1], obs[1:]):
            d0, n0, _ = cells(p, g, [y0])
            d1, n1, _ = cells(p, g, [y1])
            for o, grp in OUT_GROUPS.items():
                r = decompose(sum(d0[x] for x in grp)[0], n0[0], sum(d1[x] for x in grp)[0], n1[0])
                cum[o] += np.array([float(r['population_size']), float(r['age_sex_structure']), float(r['rates'])])
                rows.append({'geo': g, 'year': y1, 'outcome': o, 'population_size': cum[o][0],
                             'age_sex_structure': cum[o][1], 'rates': cum[o][2], 'step_from': y0})
    c = pd.DataFrame(rows)
    # carry each country's cumulative values forward over missing years before summing
    full = []
    for (g, o), s in c.groupby(['geo', 'outcome']):
        s = s.set_index('year').reindex(range(first, last + 1)).ffill()
        s['geo'], s['outcome'] = g, o
        full.append(s.reset_index())
    c = pd.concat(full, ignore_index=True)
    mal = c[c.outcome == 'MAL'].set_index(['geo', 'year'])
    comp = c[c.outcome == 'COMP'].set_index(['geo', 'year'])
    cod = (mal.rates - comp.rates).rename('rates').reset_index().assign(outcome='CODING')
    c = pd.concat([c, cod], ignore_index=True)
    c['region'] = c.geo.map(REGION)
    agg = [c]
    for keyset, name in [(['region'], None), ([], 'Europe')]:
        a = c.groupby(keyset + ['year', 'outcome'], as_index=False)[['population_size', 'age_sex_structure', 'rates']].sum()
        a['geo'] = a.region if name is None else name
        agg.append(a)
    a = c[c.geo != 'UK'].groupby(['year', 'outcome'], as_index=False)[['population_size', 'age_sex_structure', 'rates']].sum()
    agg.append(a.assign(geo='Europe excl. UK'))
    return pd.concat(agg, ignore_index=True).drop(columns=['region', 'step_from'], errors='ignore')


def set_w(p_who, grids):
    """Malignant outcome, Set W: Eurostat deaths (M, F) for 31 countries, UK from WHO."""
    e = pd.read_csv(PROCESSED / 'panel_eurostat.csv.gz')
    e = e[e.sex.isin(['M', 'F'])].assign(sex=lambda x: x.sex.map({'M': 1, 'F': 2}))
    parts = []
    for g, gr in grids.items():
        if g == 'UK':
            continue
        s = e[e.geo == g][['geo', 'year', 'sex', 'age', 'MAL', 'pop_es']].rename(columns={'pop_es': 'pop'})
        parts.append(collapse(s, gr, value_cols=['MAL', 'pop']))
    uk = p_who[p_who.geo == 'UK'][['geo', 'year', 'sex', 'age', 'MAL', 'pop']]
    w = pd.concat(parts + [uk], ignore_index=True)
    for x in ['D32', 'D33', 'D42', 'D43']:
        w[x] = 0.0
    return w


def main():
    p = load_panel()
    cs = pd.read_csv(ROOT / 'protocol' / 'country_sets.csv')
    L = cs[cs.set_L].geo.tolist()
    CC = cs[cs.set_L_complete_case].geo.tolist()
    W = cs[cs.set_W].geo.tolist()
    grids = {r.geo: [int(x) for x in r.grid.split()] for r in pd.read_csv(PROCESSED / 'audit' / 'age_grids.csv').itertuples()}
    changes = pd.read_csv(PROCESSED / 'coding_changes.csv')

    point, boot = primary_with_bootstrap(p, L)
    point.to_csv(TAB / 'decomp_primary.csv', index=False)
    rba = rate_by_age(p, L)
    chk = rba[rba.outcome == 'MAL'][['population_size', 'age_sex_structure', 'rates']].sum()
    eu = point[(point.geo == 'Europe') & (point.outcome == 'MAL')].iloc[0]
    assert np.allclose(chk.values, eu[['population_size', 'age_sex_structure', 'rates']].astype(float).values)
    rba.to_csv(TAB / 'decomp_primary_by_age_sex.csv', index=False)
    boot.to_csv(TAB / 'decomp_primary_bootstrap.csv', index=False)

    sens = []
    sens.append(run_spec(p, L, BASE_L, END_L, 'S1 narrow composite',
                         outcome_map={'MAL': ['MAL'], 'COMP': ['MAL', 'D33', 'D43'], 'COMP_NARROW': ['MAL', 'D33', 'D43'],
                                      'BU': ['D33', 'D43']}))
    p85 = pd.concat([collapse(p[p.geo == g], [0, 15] + list(range(20, 90, 5)) if g == 'NO' else list(range(0, 90, 5)),
                              value_cols=ELEM + ['pop', 'pop_who']) for g in p.geo.unique()])
    sens.append(run_spec(p85, L, BASE_L, END_L, 'S2 ages to 85+'))
    sens.append(run_spec(p, L, BASE_L, (2017, 2019), 'S3 end period 2017-2019'))
    # S7: WHO population where available for all cells of both periods
    pw = p.copy()
    ok = []
    for g in L:
        s = p[(p.geo == g) & p.year.isin(years(BASE_L) + years(END_L))]
        if g != 'UK' and s.pop_who.notna().all() and (s.pop_who > 0).all() and s.year.nunique() == 6:
            ok.append(g)
    pw['pop'] = np.where(pw.geo.isin(ok), pw.pop_who, pw['pop'])
    sens.append(run_spec(pw, L, BASE_L, END_L, 'S7 WHO population').assign(note=f'WHO population used for {len(ok)} countries: {" ".join(ok)}'))
    # S8: exclude countries with a documented or data-identified change after the first base year
    chg = changes[(changes.year > BASE_L[0]) & (changes.year <= END_L[1])]
    undated = set(pd.read_csv(PROCESSED / 'coding_changes_undated.csv').geo)      # Germany (deviations table, entry 13)
    keep8 = [g for g in L if g not in set(chg.geo) | undated]
    sens.append(run_spec(p, keep8, BASE_L, END_L, 'S8 no coding change').assign(note='countries: ' + ' '.join(keep8)))
    sens.append(run_spec(p, CC, BASE_L, END_L, 'S9 complete case'))
    # Deviations table, entry 2: Germany 2021-2022 C70-C72 from Eurostat age-sex cells
    e = pd.read_csv(PROCESSED / 'panel_eurostat.csv.gz')
    e = e[(e.geo == 'DE') & e.year.isin([2021, 2022]) & e.sex.isin(['M', 'F'])].assign(sex=lambda x: x.sex.map({'M': 1, 'F': 2}))
    e = collapse(e[['geo', 'year', 'sex', 'age', 'MAL']], grids['DE'], value_cols=['MAL'])
    # Eurostat publishes ages 0-24 as one group (stored at age 0 in panel_eurostat); it is distributed over the
    # WHO age cells 0-4 to 20-24 in proportion to the WHO deaths in those cells
    young = e.age == 0
    e_y = e[young].set_index(['geo', 'year', 'sex']).MAL
    e = e[~young]
    pde = p.merge(e.rename(columns={'MAL': 'MAL_es'}), on=['geo', 'year', 'sex', 'age'], how='left')
    m = (pde.geo == 'DE') & pde.year.isin([2021, 2022]) & (pde.age < 25)
    w_ = pde[m].groupby(['geo', 'year', 'sex']).MAL.transform('sum')
    tot = pd.Series([e_y.get((g, y, sx), np.nan) for g, y, sx in zip(pde[m].geo, pde[m].year, pde[m].sex)], index=pde[m].index)
    n_ = pde[m].groupby(['geo', 'year', 'sex']).MAL.transform('size')
    pde.loc[m, 'MAL_es'] = np.where(w_ > 0, tot * pde.loc[m, 'MAL'] / w_, tot / n_)
    pde['MAL'] = np.where(pde.MAL_es.notna(), pde.MAL_es, pde.MAL)
    sens.append(run_spec(pde, L, BASE_L, END_L, 'Germany 2021-2022 malignant from Eurostat'))
    pd.concat(sens, ignore_index=True).to_csv(TAB / 'decomp_sensitivity.csv', index=False)

    ch = chained(p, L)
    ch.to_csv(TAB / 'decomp_chained.csv', index=False)

    w = set_w(p, grids)
    rng = np.random.default_rng(SEED + 1)
    saved = dict(OUT_GROUPS)
    OUT_GROUPS.clear(); OUT_GROUPS['MAL'] = ['MAL']
    pt, bt = {}, {}
    for g in W:
        pt[g], bt[g] = country_decomp(w, g, BASE_W, END_W, rng=rng, nboot=N_BOOT)
    rows = []
    for g in W:
        rows.append({'spec': 'Set W', 'geo': g, 'outcome': 'MAL', **pt[g]['MAL']})
    dfw = pd.DataFrame(rows)
    wt = pd.concat([dfw, aggregate(dfw, True), aggregate(dfw, False)], ignore_index=True)
    brow = []
    units = {g: [g] for g in W}
    for reg in sorted(set(REGION[g] for g in W)):
        units[reg] = [g for g in W if REGION[g] == reg]
    units['Europe'] = W
    units['Europe excl. UK'] = [g for g in W if g != 'UK']
    for u, m in units.items():
        for k in ['population_size', 'age_sex_structure', 'rates']:
            v = sum(bt[g]['MAL'][k] for g in m)
            lo, hi = np.percentile(v, [2.5, 97.5])
            brow.append({'geo': u, 'outcome': 'MAL', 'component': k, 'boot_lo': lo, 'boot_hi': hi})
    # S6: Set W with WHO counts where WHO covers both periods
    whoW = []
    for g in W:
        s = p[(p.geo == g) & p.year.isin(years(BASE_W) + years(END_W))]
        if s.year.nunique() == 6 and s['pop'].notna().all():
            whoW.append(g)
    s6 = {g: country_decomp(p, g, BASE_W, END_W)[0] for g in whoW}
    s6e = {g: pt[g] for g in whoW}
    OUT_GROUPS.clear(); OUT_GROUPS.update(saved)
    s6df = pd.concat([pd.DataFrame([{'spec': 'S6 Set W WHO counts', 'geo': g, 'outcome': 'MAL', **s6[g]['MAL']} for g in whoW]),
                      pd.DataFrame([{'spec': 'S6 Set W Eurostat counts, same countries', 'geo': g, 'outcome': 'MAL', **s6e[g]['MAL']} for g in whoW])])
    s6df = pd.concat([s6df, aggregate(s6df, True)], ignore_index=True)
    wt.to_csv(TAB / 'decomp_setW.csv', index=False)
    pd.DataFrame(brow).to_csv(TAB / 'decomp_setW_bootstrap.csv', index=False)
    s6df.to_csv(TAB / 'decomp_S6_who_vs_eurostat.csv', index=False)
    eu = point[(point.geo.isin(['Europe', 'Europe excl. UK']))]
    print(eu.round(1).to_string(index=False))
    print(wt[wt.geo.isin(['Europe', 'Europe excl. UK'])].round(1).to_string(index=False))


if __name__ == '__main__':
    main()
