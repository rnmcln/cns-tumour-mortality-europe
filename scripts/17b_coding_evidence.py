"""Step 17b: analyses of the non-malignant rate component (additional analyses).

None of these analyses was specified in protocol v1.0. They test the interpretation of the 'coding
component' (rate component of the malignant outcome minus that of the composite outcome, which equals
minus the rate component of D32-D33 and D42-D43) as a transfer of deaths to malignant codes.

1. Subcodes: decomposition of D43 into D43.2 (brain, unspecified), D43.0-D43.1 (brain, supratentorial
   and infratentorial) and other D43, and of C70-C72 into C71.9 (brain, unspecified) and other malignant
   codes, in the Set L countries coding at four characters in every year of both periods.
2. Bounds: the non-malignant rate component is the largest transfer compatible with the data. A
   matched bound counts, in each country x broad age group cell, only the part of the non-malignant
   decline that is accompanied by a malignant rate increase: sum of min(max(-BU, 0), max(MAL, 0)).
   A regression estimate uses the weighted cross-country slope of the malignant on the non-malignant
   rate component (both as % of baseline composite deaths; weights = baseline composite deaths).
3. Baseline period: 2004-2006, 2005-2007 and 2006-2008 instead of 2003-2005.
4. Placebo: rate components of non-CNS D00-D48 codes and of non-CNS malignant codes (C00-C97 except
   C70-C72) in the same countries and periods.
5. Ill-defined causes: deaths at ages 75+ assigned to R00-R99 redistributed in proportion to the
   distribution of non-R deaths in each country, year, sex and age group.
6. Broad composite: composite outcome plus sellar codes (C75.1-C75.3, D35.2-D35.4, D44.3-D44.5).
7. Year-block bootstrap: years resampled with replacement within each period, then Poisson draws.
8. Country-level indicators: change in the unspecified-subcode share and in the BU/MAL ratio.
Outputs: output/tables/coding_evidence_*.csv
"""
import numpy as np
import pandas as pd

from config import OUTPUT, ROOT, SEED, BASE_L, END_L
from lib.data import load_panel
from lib.decomp import decompose

TAB = OUTPUT / 'tables'
BU = ['D32', 'D33', 'D42', 'D43']
BROAD = [(-1, 39, '0-39'), (40, 64, '40-64'), (65, 74, '65-74'), (75, 84, '75-84'), (85, 200, '85+')]


def yrs(p):
    return list(range(p[0], p[1] + 1))


def arrays(p, geo, years, cols):
    s = p[(p.geo == geo) & p.year.isin(years)].sort_values(['year', 'sex', 'age'])
    ny = len(years)
    assert s.year.nunique() == ny, (geo, years)
    nc = len(s) // ny
    d = {c: s[c].to_numpy(float).reshape(ny, nc) for c in cols}
    n = s['pop'].to_numpy(float).reshape(ny, nc)
    key = s[['sex', 'age']].iloc[:nc].to_numpy()
    return d, n, key


def comp(d0, n0, d1, n1, by_cell=False):
    return decompose(d0.mean(0) if d0.ndim == 2 else d0, n0.mean(0), d1.mean(0) if d1.ndim == 2 else d1, n1.mean(0), by_cell)


def main():
    p = load_panel()
    p = p[p['pop'].notna()].copy()
    p['BU'] = p[BU].sum(axis=1)
    p['COMP'] = p.MAL + p.BU
    p['MAL_SPEC'] = p.MAL - p.C71_9
    p['DNONCNS'] = p.DALL - p.BU
    p['CNONCNS'] = p.CALL - p.MAL
    p['COMP_BROAD'] = p.COMP + p.SELLAR
    cs = pd.read_csv(ROOT / 'protocol' / 'country_sets.csv')
    L = cs[cs.set_L].geo.tolist()
    REP = set(cs[cs.country_level_reporting].geo)
    out = {}

    # ---- 1. subcodes (four-character countries) --------------------------------------------------
    per = p[p.year.isin(yrs(BASE_L) + yrs(END_L)) & p.geo.isin(L)]
    c4 = [g for g in L if per[per.geo == g].groupby('year').MAL_3CHAR.sum().eq(0).all()
          and per[per.geo == g].icd_list.astype(str).eq('104').all()]
    cols = ['MAL', 'C71_9', 'MAL_SPEC', 'D43', 'D43_2', 'D43_01', 'D43_OTH', 'D32', 'D33', 'D42', 'BU', 'COMP']
    rows = []
    for g in c4:
        d0, n0, _ = arrays(p, g, yrs(BASE_L), cols)
        d1, n1, _ = arrays(p, g, yrs(END_L), cols)
        for c in cols:
            r = comp(d0[c], n0, d1[c], n1)
            rows.append({'geo': g, 'outcome': c, **{k: float(v) for k, v in r.items()}})
    sc = pd.DataFrame(rows)
    tot = sc.groupby('outcome', as_index=False).sum(numeric_only=True).assign(geo='C4 total')
    sc = pd.concat([sc, tot], ignore_index=True)
    sc['n_countries'] = len(c4)
    sc.to_csv(TAB / 'coding_evidence_subcodes.csv', index=False)
    out['subcodes'] = tot.set_index('outcome').rates.round(1).to_dict()

    # ---- 2. bounds -----------------------------------------------------------------------------
    cell_rows = []
    for g in L:
        d0, n0, key = arrays(p, g, yrs(BASE_L), ['MAL', 'BU', 'COMP'])
        d1, n1, _ = arrays(p, g, yrs(END_L), ['MAL', 'BU', 'COMP'])
        r = {o: comp(d0[o], n0, d1[o], n1, by_cell=True) for o in ['MAL', 'BU', 'COMP']}
        for i, (sx, age) in enumerate(key):
            b = next(lab for lo, hi, lab in BROAD if lo <= age <= hi)
            cell_rows.append({'geo': g, 'broad': b, 'MAL': r['MAL']['rates'][i], 'BU': r['BU']['rates'][i],
                              'COMP0': r['COMP']['deaths_0'][i]})
    ce = pd.DataFrame(cell_rows)
    cb = ce.groupby(['geo', 'broad'], as_index=False)[['MAL', 'BU', 'COMP0']].sum()
    cb['matched'] = np.minimum(np.clip(-cb.BU, 0, None), np.clip(cb.MAL, 0, None))
    cb.to_csv(TAB / 'coding_evidence_bounds_cells.csv', index=False)
    cc = ce.groupby('geo', as_index=False)[['MAL', 'BU', 'COMP0']].sum()
    cc['matched'] = np.minimum(np.clip(-cc.BU, 0, None), np.clip(cc.MAL, 0, None))
    cc['mal_pct'] = 100 * cc.MAL / cc.COMP0
    cc['bu_pct'] = 100 * cc.BU / cc.COMP0
    cc['reported'] = cc.geo.isin(REP)
    cc.to_csv(TAB / 'coding_evidence_country_components.csv', index=False)
    # weighted least squares slope (reported countries), bootstrap over countries for the interval
    r = cc[cc.reported].reset_index(drop=True)

    def wls(df):
        w = df.COMP0.to_numpy()
        x, y = df.bu_pct.to_numpy(), df.mal_pct.to_numpy()
        xm, ym = np.average(x, weights=w), np.average(y, weights=w)
        return np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2)
    slope = wls(r)
    rng = np.random.default_rng(SEED)
    bs = np.array([wls(r.iloc[rng.integers(0, len(r), len(r))]) for _ in range(2000)])
    bu_total = cc.BU.sum()
    from scipy.stats import spearmanr
    rho, pval = spearmanr(r.bu_pct, r.mal_pct)
    bounds = pd.DataFrame([
        {'estimate': 'upper: non-malignant rate component (coding component)', 'deaths': -bu_total},
        {'estimate': 'matched, country x broad age group', 'deaths': cb.matched.sum()},
        {'estimate': 'matched, country', 'deaths': cc.matched.sum()},
        {'estimate': 'regression: -slope x non-malignant rate component', 'deaths': slope * bu_total,
         'lo': np.percentile(bs, 2.5) * bu_total, 'hi': np.percentile(bs, 97.5) * bu_total},
    ])
    bounds['slope'] = [np.nan, np.nan, np.nan, slope]
    bounds['slope_lo'] = [np.nan, np.nan, np.nan, np.percentile(bs, 2.5)]
    bounds['slope_hi'] = [np.nan, np.nan, np.nan, np.percentile(bs, 97.5)]
    bounds['spearman_rho'] = [np.nan, np.nan, np.nan, rho]
    bounds['spearman_p'] = [np.nan, np.nan, np.nan, pval]
    bounds['n_countries'] = [len(L), len(L), len(L), len(r)]
    bounds.to_csv(TAB / 'coding_evidence_bounds.csv', index=False)
    out['bounds'] = bounds[['estimate', 'deaths']].round(0).values.tolist()

    # ---- 3. baseline period ----------------------------------------------------------------------
    rows = []
    for b0 in [2003, 2004, 2005, 2006]:
        base = (b0, b0 + 2)
        acc = {o: {'deaths_0': 0.0, 'deaths_1': 0.0, 'population_size': 0.0, 'age_sex_structure': 0.0, 'rates': 0.0}
               for o in ['MAL', 'COMP', 'D43', 'D32']}
        for g in L:
            d0, n0, _ = arrays(p, g, yrs(base), list(acc))
            d1, n1, _ = arrays(p, g, yrs(END_L), list(acc))
            for o in acc:
                for k, v in comp(d0[o], n0, d1[o], n1).items():
                    acc[o][k] += float(v)
        for o, v in acc.items():
            rows.append({'base': f'{base[0]}-{base[1]}', 'outcome': o, **v})
        rows.append({'base': f'{base[0]}-{base[1]}', 'outcome': 'CODING', 'rates': acc['MAL']['rates'] - acc['COMP']['rates']})
    pd.DataFrame(rows).to_csv(TAB / 'coding_evidence_baseline.csv', index=False)

    # ---- 4. placebo, 6. broad composite, 5. ill-defined redistribution ------------------------------
    q = p.copy()
    old = q.age >= 75
    with np.errstate(invalid='ignore', divide='ignore'):
        f = np.where(old & (q.ALL - q.R > 0), q.ALL / (q.ALL - q.R), 1.0)
    q['MAL_R'] = q.MAL * f
    q['COMP_R'] = q.COMP * f
    rows = []
    for g in L:
        cols = ['MAL', 'COMP', 'DNONCNS', 'CNONCNS', 'COMP_BROAD', 'MAL_R', 'COMP_R', 'SELLAR']
        d0, n0, _ = arrays(q, g, yrs(BASE_L), cols)
        d1, n1, _ = arrays(q, g, yrs(END_L), cols)
        for c in cols:
            rows.append({'geo': g, 'outcome': c, **{k: float(v) for k, v in comp(d0[c], n0, d1[c], n1).items()}})
    pl = pd.DataFrame(rows)
    pl = pd.concat([pl, pl.groupby('outcome', as_index=False).sum(numeric_only=True).assign(geo='Set L')], ignore_index=True)
    pl.to_csv(TAB / 'coding_evidence_placebo_sensitivity.csv', index=False)

    # ---- 7. year-block bootstrap ------------------------------------------------------------------
    rng = np.random.default_rng(SEED)
    nb = 2000
    acc = np.zeros((nb, 3))       # MAL rate, COMP rate, coding
    for g in L:
        d0, n0, _ = arrays(p, g, yrs(BASE_L), ['MAL', 'BU'])
        d1, n1, _ = arrays(p, g, yrs(END_L), ['MAL', 'BU'])
        i0 = rng.integers(0, 3, (nb, 3))
        i1 = rng.integers(0, 3, (nb, 3))
        res = {}
        for o, parts in [('MAL', ['MAL']), ('COMP', ['MAL', 'BU'])]:
            pass
        sim0 = {c: rng.poisson(d0[c][i0]).mean(1) for c in ['MAL', 'BU']}
        sim1 = {c: rng.poisson(d1[c][i1]).mean(1) for c in ['MAL', 'BU']}
        N0, N1 = n0[i0].mean(1), n1[i1].mean(1)
        rm = decompose(sim0['MAL'], N0, sim1['MAL'], N1)['rates']
        rc = decompose(sim0['MAL'] + sim0['BU'], N0, sim1['MAL'] + sim1['BU'], N1)['rates']
        acc += np.column_stack([rm, rc, rm - rc])
    yb = pd.DataFrame({'component': ['MAL rates', 'COMP rates (residual)', 'coding'],
                       'lo': np.percentile(acc, 2.5, axis=0), 'hi': np.percentile(acc, 97.5, axis=0),
                       'sd': acc.std(0)})
    yb.to_csv(TAB / 'coding_evidence_yearblock_bootstrap.csv', index=False)

    # ---- 8. country indicators ---------------------------------------------------------------------
    rows = []
    for g in L:
        s0 = p[(p.geo == g) & p.year.isin(yrs(BASE_L))]
        s1 = p[(p.geo == g) & p.year.isin(yrs(END_L))]
        four = s0.MAL_3CHAR.sum() == 0 and s1.MAL_3CHAR.sum() == 0
        rows.append({'geo': g, 'four_char_both': four,
                     'unspec_0': s0.C71_9.sum() / s0.MAL.sum() if four else np.nan,
                     'unspec_1': s1.C71_9.sum() / s1.MAL.sum() if four else np.nan,
                     'd43_2_share_0': s0.D43_2.sum() / s0.D43.sum() if four else np.nan,
                     'ratio_0': s0.BU.sum() / s0.MAL.sum(), 'ratio_1': s1.BU.sum() / s1.MAL.sum(),
                     'd43_0': s0.D43.sum() / 3, 'd43_1': s1.D43.sum() / 3, 'd32_0': s0.D32.sum() / 3, 'd32_1': s1.D32.sum() / 3})
    pd.DataFrame(rows).to_csv(TAB / 'coding_evidence_country_indicators.csv', index=False)
    print('four-character countries:', len(c4), c4)
    print(sc[sc.geo == 'C4 total'][['outcome', 'deaths_0', 'deaths_1', 'rates']].round(0).to_string(index=False))
    print(bounds.round(2).to_string(index=False))
    print(pd.read_csv(TAB / 'coding_evidence_baseline.csv').query("outcome in ['CODING','MAL','COMP']")[['base', 'outcome', 'rates']].round(0).to_string(index=False))
    print(pl[pl.geo == 'Set L'][['outcome', 'deaths_0', 'deaths_1', 'rates']].round(0).to_string(index=False))
    print(yb.round(0).to_string(index=False))


if __name__ == '__main__':
    main()
