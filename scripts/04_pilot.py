"""Step 4 (pilot): data audit, trends, decomposition, projections and backtest.

Outputs: pilot/tables/*.csv, pilot/pilot_numbers.json, pilot/figures/* (via 05_pilot_figures.py).
All results are provisional and serve to test feasibility and design choices.
"""
import json

import numpy as np
import pandas as pd

from config import PROCESSED, RAW, PILOT, COUNTRIES, ESP2013, SEED
from lib.methods import asmr, decompose_deaths, attenuated_trend_projection

T = PILOT / 'tables'
T.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(SEED)
NUM = {}

D = pd.read_csv(PROCESSED / 'deaths_long.csv', dtype={'frmat': str})
P = pd.read_csv(PROCESSED / 'population_long.csv')
A21 = list(range(0, 100, 5))                     # Eurostat groups 0-4 ... 95+
A18 = list(range(0, 90, 5))                      # harmonised groups 0-4 ... 85+
ESP21 = ESP2013[:-1] + [1500, 800, 200]
MID21 = [2.5 if a == 0 else a + 2.5 for a in A21[:-1]] + [97.0]


def collapse85(df, val):
    df = df.copy()
    df['age'] = df.age.clip(upper=85)
    keys = [c for c in df.columns if c not in (val,)]
    return df.groupby(keys, as_index=False)[val].sum()


# ---------------------------------------------------------------- 1. data audit
eu = D[D.source == 'eurostat']
who = D[D.source == 'who']
cov = []
for src, d in (('eurostat', eu), ('who', who)):
    b = d[(d.cause == 'C70-C72') & d.sex.isin(['T'] if src == 'eurostat' else ['M', 'F'])].groupby(['geo', 'year']).deaths.sum().reset_index()
    for g in COUNTRIES:
        ys = sorted(b[b.geo == g].year)
        cov.append({'source': src, 'geo': g, 'first': ys[0] if ys else None, 'last': ys[-1] if ys else None,
                    'n_years': len(ys),
                    'missing_within': ','.join(str(y) for y in range(ys[0], ys[-1] + 1) if y not in ys) if ys else ''})
cov = pd.DataFrame(cov)
cov.to_csv(T / 'audit_coverage.csv', index=False)

# WHO versus Eurostat, overlapping years (all ages, both sexes, C70-C72)
e_tot = eu[(eu.cause == 'C70-C72') & (eu.sex == 'T')].groupby(['geo', 'year']).deaths.sum()
w_tot = who[who.cause == 'C70-C72'].groupby(['geo', 'year']).deaths.sum()
agree = pd.concat([e_tot.rename('eurostat'), w_tot.rename('who')], axis=1).dropna().reset_index()
agree['diff_pct'] = (agree.who / agree.eurostat - 1) * 100
agree.to_csv(T / 'audit_who_vs_eurostat.csv', index=False)
NUM['agree_cells'] = int(len(agree))
NUM['agree_within_1pct'] = int((agree.diff_pct.abs() <= 1).sum())
NUM['agree_within_5pct'] = int((agree.diff_pct.abs() <= 5).sum())
worst = agree.assign(a=agree.diff_pct.abs()).sort_values('a', ascending=False).head(10)
NUM['agree_worst'] = worst[['geo', 'year', 'eurostat', 'who', 'diff_pct']].round(1).to_dict('records')

# certification indicators (Eurostat): ill-defined share of all deaths at 85+, by year (31 countries)
cert = eu[eu.cause.isin(['ALL', 'R00-R99', 'C70-C72', 'D00-D48']) & (eu.sex == 'T') & eu.year.between(2011, 2023)].copy()
cert['band'] = np.where(cert.age >= 85, '85+', np.where(cert.age >= 65, '65-84', '<65'))
c = cert.groupby(['year', 'band', 'cause']).deaths.sum().unstack('cause')
c['R_share_pct'] = c['R00-R99'] / c['ALL'] * 100
c['D00D48_per_C70C72'] = c['D00-D48'] / c['C70-C72']
c.round(3).to_csv(T / 'audit_certification_eurostat.csv')
# WHO: benign / uncertain CNS tumours (D32-D33, D42-D43) relative to C70-C72, by country and period
wd = who[who.cause.isin(['C70-C72', 'D32-D33_D42-D43'])].copy()
wd['period'] = pd.cut(wd.year, [1993, 2004, 2012, 2019, 2024], labels=['1994-2004', '2005-2012', '2013-2019', '2020-2024'])
r = wd.groupby(['geo', 'period', 'cause'], observed=True).deaths.sum().unstack('cause')
r['ratio_D_to_C'] = r['D32-D33_D42-D43'] / r['C70-C72']
r.round(3).to_csv(T / 'audit_benign_uncertain_ratio_who.csv')

# ------------------------------------------------- 2. trends 2011-2023 (Eurostat)
def panel(geo_list, cause='C70-C72', y0=2011, y1=2023, sexes=('T',)):
    d = eu[(eu.cause == cause) & eu.geo.isin(geo_list) & eu.sex.isin(sexes) & eu.year.between(y0, y1)]
    p = P[(P.source == 'eurostat') & P.geo.isin(geo_list) & P.sex.isin(sexes) & P.year.between(y0, y1)]
    dd = d.groupby(['year', 'age']).deaths.sum().unstack('age').reindex(columns=A21).fillna(0)
    pp = p.groupby(['year', 'age']).population.sum().unstack('age').reindex(columns=A21)
    return dd, pp


def complete_geos(y0, y1, cause='C70-C72'):
    ok = []
    for g in COUNTRIES:
        ys = set(eu[(eu.geo == g) & (eu.cause == cause)].year)
        ps = set(P[(P.source == 'eurostat') & (P.geo == g)].year)
        if all(y in ys and y in ps for y in range(y0, y1 + 1)):
            ok.append(g)
    return ok


GEOS = complete_geos(2011, 2023)
NUM['n_countries_complete_2011_2023'] = len(GEOS)
NUM['countries_incomplete_2011_2023'] = sorted(set(COUNTRIES) - set(GEOS))

trend_rows = []
for label, geos in [('Europe-' + str(len(GEOS)), GEOS)] + [(g, [g]) for g in GEOS]:
    for sexes, sl in ((('T',), 'both'), (('M',), 'male'), (('F',), 'female')):
        dd, pp = panel(geos, sexes=sexes)
        for y in dd.index:
            r_, se = asmr(dd.loc[y].values, pp.loc[y].values, ESP21)
            trend_rows.append({'unit': label, 'sex': sl, 'year': y, 'deaths': dd.loc[y].sum(),
                               'population': pp.loc[y].sum(), 'asmr': r_, 'asmr_se': se,
                               'deaths_75plus': dd.loc[y, [a for a in A21 if a >= 75]].sum()})
trend = pd.DataFrame(trend_rows)
trend.round(4).to_csv(T / 'trends_eurostat.csv', index=False)


def apc(unit, sex, y0, y1):
    """Annual percentage change in ASMR: weighted log-linear regression of ASMR on year."""
    s = trend[(trend.unit == unit) & (trend.sex == sex) & trend.year.between(y0, y1)]
    w = (s.asmr / s.asmr_se) ** 2
    X = np.column_stack([np.ones(len(s)), s.year - y0])
    y = np.log(s.asmr)
    W = np.diag(w)
    b = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
    res = y - X @ b
    phi = max(1.0, float((w * res ** 2).sum() / (len(s) - 2)))
    cov_ = np.linalg.inv(X.T @ W @ X) * phi
    se = np.sqrt(cov_[1, 1])
    return 100 * (np.exp(b[1]) - 1), 100 * (np.exp(b[1] - 1.96 * se) - 1), 100 * (np.exp(b[1] + 1.96 * se) - 1)


apc_rows = []
for unit in trend.unit.unique():
    for sex in ('both', 'male', 'female'):
        for (y0, y1) in ((2011, 2019), (2011, 2023)):
            e, lo, hi = apc(unit, sex, y0, y1)
            apc_rows.append({'unit': unit, 'sex': sex, 'period': f'{y0}-{y1}', 'apc': e, 'lo': lo, 'hi': hi})
apc_df = pd.DataFrame(apc_rows)
apc_df.round(2).to_csv(T / 'trends_apc_eurostat.csv', index=False)
EU = 'Europe-' + str(len(GEOS))
s = trend[(trend.unit == EU) & (trend.sex == 'both')].set_index('year')
NUM['eu_deaths_2011'] = int(s.loc[2011, 'deaths']); NUM['eu_deaths_2023'] = int(s.loc[2023, 'deaths'])
NUM['eu_asmr_2011'] = round(s.loc[2011, 'asmr'], 2); NUM['eu_asmr_2023'] = round(s.loc[2023, 'asmr'], 2)
NUM['eu_share_75plus_2011'] = round(s.loc[2011, 'deaths_75plus'] / s.loc[2011, 'deaths'] * 100, 1)
NUM['eu_share_75plus_2023'] = round(s.loc[2023, 'deaths_75plus'] / s.loc[2023, 'deaths'] * 100, 1)
NUM['eu_apc'] = apc_df[apc_df.unit == EU].round(2).to_dict('records')

# age-specific rates, Europe aggregate, 2011-2013 vs 2021-2023
dd, pp = panel(GEOS)
r0 = dd.loc[2011:2013].sum() / pp.loc[2011:2013].sum() * 1e5
r1 = dd.loc[2021:2023].sum() / pp.loc[2021:2023].sum() * 1e5
pd.DataFrame({'age': A21, 'rate_2011_2013': r0.values, 'rate_2021_2023': r1.values,
              'pct_change': (r1 / r0 - 1).values * 100}).round(3).to_csv(T / 'age_specific_rates_eu.csv', index=False)

# long series (WHO), countries with complete 2000-2022 C70-C72 and population
who_b = who[who.cause == 'C70-C72']
wp = P[P.source == 'who']
long_geos = []
for g in COUNTRIES:
    ys = set(who_b[who_b.geo == g].year); ps = set(wp[wp.geo == g].year)
    if all(y in ys and y in ps for y in range(2000, 2023)):
        long_geos.append(g)
NUM['who_complete_2000_2022'] = long_geos
lrows = []
for g in long_geos:
    d = who_b[(who_b.geo == g) & who_b.year.between(2000, 2022)].groupby(['year', 'age']).deaths.sum().unstack().reindex(columns=A18).fillna(0)
    p = wp[(wp.geo == g) & wp.year.between(2000, 2022)].groupby(['year', 'age']).population.sum().unstack().reindex(columns=A18)
    for y in d.index:
        r_, se = asmr(d.loc[y].values, p.loc[y].values, ESP2013)
        lrows.append({'geo': g, 'year': y, 'deaths': d.loc[y].sum(), 'asmr': r_, 'asmr_se': se})
pd.DataFrame(lrows).round(4).to_csv(T / 'trends_who_long.csv', index=False)

# ----------------------------------------------- 3. reference life table for YLL
le = pd.read_csv(RAW / 'eurostat' / 'demo_mlexpec.csv.gz').rename(columns={'geo\\TIME_PERIOD': 'geo'})
ch = le[(le.geo == 'CH') & (le.sex == 'T')][['age', '2023']].dropna()


def a_num(a):
    return 0 if a == 'Y_LT1' else (int(a[1:]) if a[1:].isdigit() else (95 if a == 'Y_GE95' else None))
ch['a'] = ch.age.map(a_num)
ch = ch.dropna(subset=['a']).sort_values('a')
REF_LE = np.interp(MID21, ch.a.values, ch['2023'].values)
NUM['reference_life_table'] = f'Switzerland 2023, both sexes; ages {int(ch.a.min())}-{int(ch.a.max())}'
pd.DataFrame({'age': A21, 'midpoint': MID21, 'remaining_le': REF_LE}).round(3).to_csv(T / 'reference_le.csv', index=False)

# ------------------------------------------------------- 4. decomposition 2011-2013 to 2021-2023
dec_rows = []
for label, geos in [(EU, GEOS)] + [(g, [g]) for g in GEOS]:
    dd, pp = panel(geos)
    D0, D1 = dd.loc[2011:2013].mean(), dd.loc[2021:2023].mean()
    N0, N1 = pp.loc[2011:2013].mean(), pp.loc[2021:2023].mean()
    for measure, w in (('deaths', None), ('yll', REF_LE)):
        res = decompose_deaths(N0.values, (D0 / N0).values, N1.values, (D1 / N1).values, weight=w)
        # parametric bootstrap of death counts (Poisson), 1000 replicates
        boots = []
        for _ in range(1000):
            d0 = rng.poisson(dd.loc[2011:2013].values).mean(axis=0)
            d1 = rng.poisson(dd.loc[2021:2023].values).mean(axis=0)
            rb = decompose_deaths(N0.values, d0 / N0.values, N1.values, d1 / N1.values, weight=w)
            boots.append([rb['population_size'], rb['age_structure'], rb['age_specific_rates']])
        q = np.percentile(np.array(boots), [2.5, 97.5], axis=0)
        dec_rows.append({'unit': label, 'measure': measure, 'base_2011_2013': res['total_0'],
                         'latest_2021_2023': res['total_1'], 'change': res['total_1'] - res['total_0'],
                         'population_size': res['population_size'], 'age_structure': res['age_structure'],
                         'age_specific_rates': res['age_specific_rates'],
                         'rates_lo': q[0, 2], 'rates_hi': q[1, 2]})
dec = pd.DataFrame(dec_rows)
dec.round(1).to_csv(T / 'decomposition_2011_13_to_2021_23.csv', index=False)
NUM['decomposition_eu'] = dec[dec.unit == EU].round(1).to_dict('records')

# exploratory: do country differences in rate trends track shifts between malignant and
# benign/uncertain CNS codes? (WHO D/C ratio 2011-2013 vs 2020-2022, same years as decomposition)
wr = wd[wd.year.between(2011, 2013) | wd.year.between(2020, 2022)].copy()
wr['p'] = np.where(wr.year <= 2013, 'p0', 'p1')
wr = wr.groupby(['geo', 'p', 'cause']).deaths.sum().unstack('cause')
wr['ratio'] = wr['D32-D33_D42-D43'] / wr['C70-C72']
wr = wr.ratio.unstack('p').dropna()
x = dec[(dec.measure == 'deaths') & (dec.unit != EU)].set_index('unit')
x = x.assign(rates_pct=x.age_specific_rates / x.base_2011_2013 * 100).join(wr, how='inner')
x['ratio_change'] = x.p1 - x.p0
x = x[x.base_2011_2013 >= 100]                    # exclude countries with <100 deaths per year
from scipy.stats import spearmanr
rho, pval = spearmanr(x.ratio_change, x.rates_pct)
NUM['explore_ratio_vs_rates'] = {'n_countries': int(len(x)), 'spearman_rho': round(float(rho), 2),
                                 'p_value': round(float(pval), 3)}
x[['rates_pct', 'p0', 'p1', 'ratio_change']].round(3).to_csv(T / 'explore_coding_shift_vs_rates.csv')

# sensitivity: malignant only versus malignant + benign/uncertain CNS codes (WHO deaths, Eurostat
# populations), countries with complete WHO data in 2011-2013 and 2020-2022
pe = P[(P.source == 'eurostat') & P.sex.isin(['M', 'F'])].copy()
pe['age'] = pe.age.clip(upper=85)
pe = pe.groupby(['geo', 'year', 'age'], as_index=False).population.sum()
yrs = list(range(2011, 2014)) + list(range(2020, 2023))
wc = who[who.cause.isin(['C70-C72', 'D32-D33_D42-D43']) & who.year.isin(yrs)]
ok = [g for g in COUNTRIES if all(((wc.geo == g) & (wc.year == y)).any() for y in yrs) and g != 'DE']
sens = []
for defn, causes in (('C70-C72', ['C70-C72']), ('C70-C72 + D32-D33 + D42-D43', ['C70-C72', 'D32-D33_D42-D43'])):
    d = wc[wc.geo.isin(ok) & wc.cause.isin(causes)].groupby(['year', 'age']).deaths.sum().unstack().reindex(columns=A18).fillna(0)
    p = pe[pe.geo.isin(ok) & pe.year.isin(yrs)].groupby(['year', 'age']).population.sum().unstack().reindex(columns=A18)
    N0, N1 = p.loc[2011:2013].mean(), p.loc[2020:2022].mean()
    D0, D1 = d.loc[2011:2013].mean(), d.loc[2020:2022].mean()
    res = decompose_deaths(N0.values, (D0 / N0).values, N1.values, (D1 / N1).values)
    sens.append({'definition': defn, 'n_countries': len(ok), **{k: round(v, 1) for k, v in res.items()},
                 'rates_pct_of_baseline': round(res['age_specific_rates'] / res['total_0'] * 100, 2)})
pd.DataFrame(sens).to_csv(T / 'sensitivity_case_definition_who.csv', index=False)
NUM['sensitivity_case_definition'] = sens
NUM['sensitivity_case_definition_countries'] = ok

# ---------------------------------------------------------- 5. projections to 2040
proj = pd.read_csv(RAW / 'eurostat' / 'proj_23np.csv.gz').rename(columns={'geo\\TIME_PERIOD': 'geo'})
proj = proj[proj.sex.isin(['M', 'F']) & proj.geo.isin(GEOS)]
pyears = [c for c in proj.columns if c.isdigit() and 2024 <= int(c) <= 2040]


def single_age(a):
    if a == 'Y_LT1':
        return 0
    if a == 'Y_GE100':
        return 100
    return int(a[1:]) if a.startswith('Y') and a[1:].isdigit() else None
proj['a1'] = proj.age.map(single_age)
proj = proj[proj.a1.notna()]
proj['age'] = (proj.a1 // 5 * 5).clip(upper=95).astype(int)
PJ = proj.melt(id_vars=['geo', 'sex', 'age', 'projection'], value_vars=pyears, var_name='year', value_name='population')
PJ['year'] = PJ.year.astype(int)
PJ = PJ.groupby(['projection', 'geo', 'year', 'age'], as_index=False).population.sum()

PGEOS = sorted(set(GEOS) & set(PJ.geo.unique()))   # EUROPOP2023 has no projection for Serbia
NUM['projection_countries'] = len(PGEOS)
PJ = PJ[PJ.geo.isin(PGEOS)]
dd, pp = panel(PGEOS)
proj_rows = []
for variant in sorted(PJ.projection.unique()):
    pv = PJ[PJ.projection == variant].groupby(['year', 'age']).population.sum().unstack().reindex(columns=A21)
    for base_label, (b0, b1) in (('rates 2017-2019', (2017, 2019)), ('rates 2021-2023', (2021, 2023))):
        rate = dd.loc[b0:b1].sum() / pp.loc[b0:b1].sum()
        for y in pv.index:
            dy = (pv.loc[y] * rate)
            proj_rows.append({'variant': variant, 'scenario': f'A: constant {base_label}', 'year': y,
                              'deaths': dy.sum(), 'yll': (dy * REF_LE).sum()})
    if variant == 'BSL':
        for fit_from, lab in ((2011, 'B: attenuated trend fitted 2011-2023'), (2011, None)):
            if lab is None:
                continue
            out, slopes = attenuated_trend_projection(dd.index.values, dd.values, pp.values,
                                                      pv.index.values, pv.values, fit_from=fit_from)
            NUM['trend_slopes_pct_per_year'] = dict(zip([str(a) for a in A21], np.round((np.exp(slopes) - 1) * 100, 2)))
            for i, y in enumerate(pv.index):
                proj_rows.append({'variant': variant, 'scenario': lab, 'year': y,
                                  'deaths': out[i].sum(), 'yll': (out[i] * REF_LE).sum()})
for y in dd.index:
    proj_rows.append({'variant': 'observed', 'scenario': 'observed', 'year': y,
                      'deaths': dd.loc[y].sum(), 'yll': (dd.loc[y].values * REF_LE).sum()})
pr = pd.DataFrame(proj_rows)
pr.round(1).to_csv(T / 'projections_eu.csv', index=False)
obs23 = dd.loc[2021:2023].sum(axis=1).mean()
NUM['proj_countries_deaths_mean_2021_2023'] = round(float(obs23), 1)
NUM['proj_countries_yll_mean_2021_2023'] = round(float((dd.loc[2021:2023].values * REF_LE).sum(axis=1).mean()), 1)
NUM['projection_2040'] = pr[pr.year == 2040].round(0).to_dict('records')

# --------------------------------------------------------------------- 6. backtest
dd, pp = panel(GEOS)
bt = []
for fit_to in (2017,):
    sel = dd.index <= fit_to
    out, _ = attenuated_trend_projection(dd.index.values[sel], dd.values[sel], pp.values[sel],
                                         dd.index.values[~sel], pp.values[~sel], fit_from=2011)
    rate_c = dd.loc[fit_to - 2:fit_to].sum() / pp.loc[fit_to - 2:fit_to].sum()
    for i, y in enumerate(dd.index.values[~sel]):
        bt.append({'fit_to': fit_to, 'year': y, 'observed': dd.loc[y].sum(),
                   'trend_model': out[i].sum(), 'constant_rates': (pp.loc[y] * rate_c).sum()})
# WHO long series backtest: fit 2000-2011, project 2012-2022 (countries with complete series)
for g in long_geos:
    d = who_b[(who_b.geo == g) & who_b.year.between(2000, 2022)].groupby(['year', 'age']).deaths.sum().unstack().reindex(columns=A18).fillna(0)
    p = wp[(wp.geo == g) & wp.year.between(2000, 2022)].groupby(['year', 'age']).population.sum().unstack().reindex(columns=A18)
    sel = d.index <= 2011
    out, _ = attenuated_trend_projection(d.index.values[sel], d.values[sel], p.values[sel],
                                         d.index.values[~sel], p.values[~sel], fit_from=2000)
    rc = d.loc[2009:2011].sum() / p.loc[2009:2011].sum()
    for i, y in enumerate(d.index.values[~sel]):
        bt.append({'fit_to': 2011, 'unit': g, 'year': y, 'observed': d.loc[y].sum(),
                   'trend_model': out[i].sum(), 'constant_rates': (p.loc[y] * rc).sum()})
# pooled-slope variant: age-specific slopes estimated on the pooled long-series countries,
# applied to each country's own fitted level (shrinks noisy country slopes to the European slope)
LD, LP = {}, {}
for g in long_geos:
    LD[g] = who_b[(who_b.geo == g) & who_b.year.between(2000, 2022)].groupby(['year', 'age']).deaths.sum().unstack().reindex(columns=A18).fillna(0)
    LP[g] = wp[(wp.geo == g) & wp.year.between(2000, 2022)].groupby(['year', 'age']).population.sum().unstack().reindex(columns=A18)
pool_d = sum(LD[g].loc[2000:2011] for g in long_geos)
pool_p = sum(LP[g].loc[2000:2011] for g in long_geos)
_, pooled_slopes = attenuated_trend_projection(pool_d.index.values, pool_d.values, pool_p.values,
                                               [2012], pool_p.values[-1:], fit_from=2000)
for g in long_geos:
    d, p = LD[g], LP[g]
    sel = d.index <= 2011
    rc = d.loc[2009:2011].sum() / p.loc[2009:2011].sum()
    for i, y in enumerate(d.index.values[~sel]):
        h = y - 2011
        eff = sum(1 - (0.0, 0.25, 0.5, 0.75)[min((k - 1) // 5, 3)] for k in range(1, h + 1))
        pooled = (p.loc[y] * rc * np.exp(pooled_slopes * eff)).sum()
        for row in bt:
            if row.get('unit') == g and row['year'] == y:
                row['pooled_slope_model'] = pooled
bt = pd.DataFrame(bt)
bt['unit'] = bt.get('unit').fillna(EU) if 'unit' in bt else EU
bt.round(1).to_csv(T / 'backtest.csv', index=False)
for m in ('trend_model', 'constant_rates', 'pooled_slope_model'):
    bt[f'ape_{m}'] = (bt[m] / bt.observed - 1).abs() * 100
NUM['backtest_eu_2018_2023_mape'] = bt[bt.unit == EU][['ape_trend_model', 'ape_constant_rates']].mean().round(1).to_dict()
NUM['backtest_who_2012_2022_mape_median_over_countries'] = (
    bt[bt.unit != EU].groupby('unit')[['ape_trend_model', 'ape_constant_rates', 'ape_pooled_slope_model']].mean().median().round(1).to_dict())
NUM['backtest_who_2012_2022_signed_error_pooled_countries'] = {
    m: round(float((bt[bt.unit != EU].groupby('year')[m].sum() / bt[bt.unit != EU].groupby('year')['observed'].sum() - 1).mean() * 100), 1)
    for m in ('trend_model', 'constant_rates', 'pooled_slope_model')}

(PILOT / 'pilot_numbers.json').write_text(json.dumps(NUM, indent=1, default=str))
print(json.dumps(NUM, indent=1, default=str)[:6000])
