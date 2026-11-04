"""Step 15: projection models, rolling-origin evaluation and projections (protocol v1.0, section 5.4).

Run in two phases around the R step 16 (BAPC):
  python 15_projections.py prepare    writes data/processed/bapc_jobs.csv.gz
  Rscript 16_bapc.R                   writes output/models/bapc_results.csv.gz
  python 15_projections.py evaluate   backtest metrics, model selection, projections

Backtest: Set L countries with WHO data from 2000 (the panel starts in 2000); origins 2009, 2011, 2013,
2015; horizons 5 and 10 years; observed person-years of the target years; outcomes MAL and COMP; models
M1 constant rates, M2 damped log-linear, M3 age-drift (lib/projmodels.py) and M5 BAPC (R-INLA).
Nordpred (model 4) was not run (deviations table, entry 3). Metrics per country, origin and horizon: absolute and
signed percentage error of total deaths and of deaths at 70+, and coverage of the 95% prediction
interval; target years 2020-2022 are reported separately and excluded from model selection.
Selection rule (protocol): a trend model replaces M1 only if its median absolute percentage error of
total deaths is lower at both horizons; among several, the one with the lower mean of the two medians.
Projections: origin 2019 for all Set L countries (M1 rates 2015-2019; trend models fitted 2010-2019 or,
for BAPC, 2000-2019), horizons 2040 (M1) and 2035 (trend models); EUROPOP2023 baseline and variants
(step 14). Malignant outcome for all 32 countries (Set W data, M1). Sensitivity S5 (rates of the latest
three years, 2021-2023 or 2020-2022 where 2023 is not available) and S11 (population variants; UN WPP
denominators for the UK).
Uncertainty is reported separately: Poisson (variance = expected count), parameter (draws of the
expected count), scenario (population variants; models).
"""
import sys

import numpy as np
import pandas as pd

from config import OUTPUT, PROCESSED, RAW, ROOT, SEED, REGION
from lib.data import load_panel, collapse
from lib.projmodels import m1_constant, m2_damped, m3_agedrift

TAB = OUTPUT / 'tables'
MOD = OUTPUT / 'models'
ORIGINS = [2009, 2011, 2013, 2015]
HORIZONS = [5, 10]
OUTCOMES = ['MAL', 'COMP']
NSIM = 1000
COVID = {2020, 2021, 2022}


def setup():
    p = load_panel()
    cs = pd.read_csv(ROOT / 'protocol' / 'country_sets.csv')
    L = cs[cs.set_L].geo.tolist()
    p = p[p.geo.isin(L) & p['pop'].notna()]
    first = p.groupby('geo').year.min()
    bt = [g for g in L if first[g] <= 2000]
    return p, cs, L, bt


def arrays(p, geo, outcome):
    s = p[p.geo == geo].sort_values(['year', 'sex', 'age'])
    yrs = np.array(sorted(s.year.unique()))
    nc = len(s) // len(yrs)
    D = s[outcome].to_numpy(float).reshape(len(yrs), nc)
    N = s['pop'].to_numpy(float).reshape(len(yrs), nc)
    cells = s[['sex', 'age']].iloc[:nc].reset_index(drop=True)
    assert (s.groupby('year').size() == nc).all()
    return yrs, D, N, cells


def prepare():
    p, cs, L, bt = setup()
    pp = pd.read_csv(PROCESSED / 'pop_projections.csv.gz')
    pp = pp[pp.variant == 'BSL'].set_index(['geo', 'year', 'sex', 'age']).py.to_dict()
    jobs = []
    for geo in L:
        for o in OUTCOMES:
            yrs, D, N, cells = arrays(p, geo, o)
            ad = cells.age >= 15
            origins = (ORIGINS if geo in bt else []) + [2019]
            for origin in origins:
                hist = np.arange(2000 if yrs.min() <= 2000 else yrs.min(), origin + 1)
                horizon = 10 if origin < 2019 else 16
                pred = np.arange(origin + 1, origin + horizon + 1)
                for s in (1, 2):
                    cm = (cells.sex == s) & ad
                    for y in np.r_[hist, pred]:
                        for c in np.where(cm)[0]:
                            age = int(cells.age[c])
                            if y <= origin:
                                if y in yrs:
                                    d, n = D[yrs == y][0, c], N[yrs == y][0, c]
                                else:                     # missing year in the history (Slovakia)
                                    lo, hi = yrs[yrs < y].max(), yrs[yrs > y].min()
                                    n = np.sqrt(N[yrs == lo][0, c] * N[yrs == hi][0, c])
                                    d = np.nan
                            else:
                                d = np.nan
                                if origin == 2019:
                                    n = float(pp[(geo, y, s, age)])
                                elif y in yrs:
                                    n = N[yrs == y][0, c]
                                else:
                                    n = N[yrs == yrs.max()][0, c]
                            jobs.append((geo, o, origin, s, y, age, d, n))
    j = pd.DataFrame(jobs, columns=['geo', 'outcome', 'origin', 'sex', 'year', 'age', 'deaths', 'pop'])
    j.to_csv(PROCESSED / 'bapc_jobs.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
    print('BAPC jobs:', j.groupby(['geo', 'outcome', 'origin', 'sex']).ngroups)


def run_models(D, N, yrs, origin, Nt, tyears, cells, rng, bapc=None, geo=None, outcome=None):
    """Returns {model: (mean [nt, nc], draws_expected [nsim, nt, nc] or None, sd_total info)}"""
    out = {}
    m, d = m1_constant(D, N, yrs, origin, Nt, rng, NSIM)
    out['M1 constant'] = (m, d)
    m, d, _ = m2_damped(D, N, yrs, origin, Nt, tyears, rng, NSIM)
    out['M2 damped log-linear'] = (m, d)
    m, d, _ = m3_agedrift(D, N, yrs, origin, Nt, tyears, cells.sex.to_numpy(), rng, NSIM)
    out['M3 age-drift'] = (m, d)
    if bapc is not None:
        b = bapc[(bapc.geo == geo) & (bapc.outcome == outcome) & (bapc.origin == origin)]
        if len(b):
            m1m, m1d = out['M1 constant']
            mean = m1m.copy()
            var = np.zeros_like(mean)          # BAPC predictive variance (includes Poisson noise)
            ok = True
            for c in range(len(cells)):
                if cells.age[c] < 15:
                    var[:, c] = np.var(m1d[:, :, c], axis=0) + m1m[:, c]
                    continue
                r = b[(b.sex == cells.sex[c]) & (b.age == cells.age[c])].set_index('year').reindex(tyears)
                if r['mean'].isna().all():
                    ok = False
                    break
                # years beyond the BAPC horizon (2036-2040 in the final projection) stay missing
                mean[:, c] = r['mean'].to_numpy() * (Nt[:, c] / r['pop'].to_numpy())
                var[:, c] = (r['sd'].to_numpy() * (Nt[:, c] / r['pop'].to_numpy())) ** 2
            if ok:
                out['M5 BAPC'] = (mean, var)
    return out


def summarise(mean, draws, rng, mask=None):
    """Point, Poisson-noised 95% PI and parameter-only 95% interval of the total over cells (mask)."""
    if mask is None:
        mask = np.ones(mean.shape[1], bool)
    pt = mean[:, mask].sum(1)
    if draws.ndim == 2:                         # BAPC: variance array -> normal approximation
        sd = np.sqrt(draws[:, mask].sum(1))
        par_sd = np.sqrt(np.clip(sd ** 2 - pt, 0, None))
        return pt, pt - 1.96 * sd, pt + 1.96 * sd, pt - 1.96 * par_sd, pt + 1.96 * par_sd
    e = draws[:, :, mask].sum(2)
    y = rng.poisson(e)
    return (pt, np.percentile(y, 2.5, 0), np.percentile(y, 97.5, 0), np.percentile(e, 2.5, 0), np.percentile(e, 97.5, 0))


def backtest(p, bt, bapc):
    rng = np.random.default_rng(SEED + 10)
    rows = []
    for geo in bt:
        for o in OUTCOMES:
            yrs, D, N, cells = arrays(p, geo, o)
            old = (cells.age >= 70).to_numpy()
            o85m = (cells.age >= 85).to_numpy()      # deaths at 85+ (additional analysis)
            for origin in ORIGINS:
                tyears = [origin + h for h in HORIZONS if origin + h in yrs]
                if not tyears:
                    continue
                Nt = np.vstack([N[yrs == t][0] for t in tyears])
                hist = yrs <= origin
                res = run_models(D[hist], N[hist], yrs[hist], origin, Nt, np.array(tyears), cells, rng, bapc, geo, o)
                for mname, (mean, dr) in res.items():
                    tot = summarise(mean, dr, rng)
                    o70 = summarise(mean, dr, rng, old)
                    o85 = summarise(mean, dr, rng, o85m)
                    for i, t in enumerate(tyears):
                        obs = D[yrs == t][0].sum()
                        obs70 = D[yrs == t][0][old].sum()
                        rows.append({'geo': geo, 'outcome': o, 'model': mname, 'origin': origin, 'horizon': t - origin,
                                     'target': t, 'observed': obs, 'predicted': tot[0][i], 'pi_lo': tot[1][i], 'pi_hi': tot[2][i],
                                     'observed_70': obs70, 'predicted_70': o70[0][i],
                                     'observed_85': D[yrs == t][0][o85m].sum(), 'predicted_85': o85[0][i]})
    b = pd.DataFrame(rows)
    b['pe'] = 100 * (b.predicted - b.observed) / b.observed
    b['ape'] = b.pe.abs()
    b['ape_70'] = 100 * (b.predicted_70 - b.observed_70).abs() / b.observed_70
    b['pe_70'] = 100 * (b.predicted_70 - b.observed_70) / b.observed_70
    b['pe_85'] = 100 * (b.predicted_85 - b.observed_85) / b.observed_85
    b['covered'] = (b.observed >= b.pi_lo) & (b.observed <= b.pi_hi)
    b['covid_target'] = b.target.isin(COVID)
    return b


def select(b):
    m = b.groupby(['covid_target', 'outcome', 'model', 'horizon']).agg(
        n=('ape', 'size'), median_ape=('ape', 'median'), mean_ape=('ape', 'mean'), median_pe=('pe', 'median'),
        median_ape_70=('ape_70', 'median'), median_pe_70=('pe_70', 'median'), median_pe_85=('pe_85', 'median'),
        coverage=('covered', 'mean')).reset_index()
    sel = []
    for o in OUTCOMES:
        s = m[(~m.covid_target) & (m.outcome == o)].pivot(index='model', columns='horizon', values='median_ape')
        base = s.loc['M1 constant']
        better = [mm for mm in s.index if mm != 'M1 constant' and (s.loc[mm] < base).all()]
        chosen = min(better, key=lambda mm: s.loc[mm].mean()) if better else 'M1 constant'
        sel.append({'outcome': o, 'selected': chosen, 'qualifying_trend_models': '; '.join(better) or 'none',
                    **{f'M1_median_ape_h{h}': base[h] for h in HORIZONS}})
    return m, pd.DataFrame(sel)


def europe_backtest(b):
    """Aggregate over backtest countries (sum of predicted and observed) per model, origin, horizon."""
    g = b.groupby(['outcome', 'model', 'origin', 'horizon', 'target'], as_index=False)[
        ['observed', 'predicted', 'observed_70', 'predicted_70', 'observed_85', 'predicted_85']].sum()
    n = b.groupby(['outcome', 'model', 'origin', 'horizon', 'target']).geo.nunique().rename('n_countries').reset_index()
    g = g.merge(n)
    g['pe'] = 100 * (g.predicted - g.observed) / g.observed
    g['pe_70'] = 100 * (g.predicted_70 - g.observed_70) / g.observed_70
    g['pe_85'] = 100 * (g.predicted_85 - g.observed_85) / g.observed_85
    return g


def project(p, L, bapc, pp):
    """Origin 2019 projections for Set L, both outcomes (+ BU from M1), all models; variants for M1."""
    rng = np.random.default_rng(SEED + 20)
    rows, draws_store = [], {}
    for geo in L:
        for o in OUTCOMES + ['BU']:
            yrs, D, N, cells = arrays(p, geo, o)
            hist = yrs <= 2019
            for v in sorted(pp.variant.unique()):
                q = pp[(pp.geo == geo) & (pp.variant == v)].sort_values(['year', 'sex', 'age'])
                tyears = np.array(sorted(q.year.unique()))
                Nt = q.py.to_numpy().reshape(len(tyears), -1)
                qc = q[q.year == tyears[0]][['sex', 'age']].reset_index(drop=True)
                assert (qc.values == cells.values).all(), geo
                if v == 'BSL' and o != 'BU':
                    res = run_models(D[hist], N[hist], yrs[hist], 2019, Nt, tyears, cells, rng, bapc, geo, o)
                else:
                    m, d = m1_constant(D[hist], N[hist], yrs[hist], 2019, Nt, rng, NSIM)
                    res = {'M1 constant': (m, d)}
                for mname, (mean, dr) in res.items():
                    keep = tyears <= (2040 if mname == 'M1 constant' else 2035)
                    tot = summarise(mean, dr, rng)
                    o70 = mean[:, (cells.age >= 70).to_numpy()].sum(1)
                    o85 = mean[:, (cells.age >= 85).to_numpy()].sum(1)
                    for i, t in enumerate(tyears):
                        if not keep[i]:
                            continue
                        rows.append({'geo': geo, 'outcome': o, 'model': mname, 'variant': v, 'year': t,
                                     'deaths': tot[0][i], 'pi_lo': tot[1][i], 'pi_hi': tot[2][i],
                                     'par_lo': tot[3][i], 'par_hi': tot[4][i], 'deaths_70': o70[i], 'deaths_85': o85[i]})
                    if dr.ndim == 3:
                        draws_store[(geo, o, mname, v)] = (tyears[keep], dr[:, keep].sum(2))
    pr = pd.DataFrame(rows)
    # BAPC stability flag: predictive sd of the 2035 country total larger than its mean
    b35 = pr[(pr.model == 'M5 BAPC') & (pr.year == 2035)]
    unstable = set(map(tuple, b35[(b35.pi_hi - b35.pi_lo) / 3.92 > b35.deaths][['geo', 'outcome']].values))
    pd.DataFrame(sorted(unstable), columns=['geo', 'outcome']).to_csv(TAB / 'bapc_unstable.csv', index=False)
    # Europe and regional totals with parameter uncertainty from summed draws (independent across countries)
    agg = []
    units = {'Set L': L, 'Set L excl. UK': [g for g in L if g != 'UK']}
    for r in sorted(set(REGION[g] for g in L)):
        units[r] = [g for g in L if REGION[g] == r]
    for (o, mname, v), s in pr.groupby(['outcome', 'model', 'variant']):
        uu = dict(units)
        if mname == 'M5 BAPC':
            uu['Set L, stable BAPC fits'] = [g for g in L if (g, o) not in unstable]
        for u, mem in uu.items():
            t = s[s.geo.isin(mem)].groupby('year')[['deaths', 'deaths_70', 'deaths_85']].sum()
            if mname == 'M5 BAPC':
                ss = s[s.geo.isin(mem)].copy()
                var = (((ss.pi_hi - ss.pi_lo) / 3.92) ** 2).groupby(ss.year).sum()
                t['pi_lo'], t['pi_hi'] = t.deaths - 1.96 * np.sqrt(var), t.deaths + 1.96 * np.sqrt(var)
                pv = np.clip(var - t.deaths, 0, None)
                t['par_lo'], t['par_hi'] = t.deaths - 1.96 * np.sqrt(pv), t.deaths + 1.96 * np.sqrt(pv)
            else:
                e = sum(draws_store[(g, o, mname, v)][1] for g in mem)
                yrs_ = draws_store[(mem[0], o, mname, v)][0]
                y = rng.poisson(e)
                t['par_lo'] = pd.Series(np.percentile(e, 2.5, 0), index=yrs_)
                t['par_hi'] = pd.Series(np.percentile(e, 97.5, 0), index=yrs_)
                t['pi_lo'] = pd.Series(np.percentile(y, 2.5, 0), index=yrs_)
                t['pi_hi'] = pd.Series(np.percentile(y, 97.5, 0), index=yrs_)
            t['poisson_lo'] = t.deaths - 1.96 * np.sqrt(t.deaths)
            t['poisson_hi'] = t.deaths + 1.96 * np.sqrt(t.deaths)
            agg.append(t.reset_index().assign(geo=u, outcome=o, model=mname, variant=v))
    return pr, pd.concat(agg, ignore_index=True)


def latest_rates_projection(p, L, pp):
    """S5: constant rates of the latest three observed years (2021-2023, or 2020-2022)."""
    rows = []
    for geo in L:
        for o in OUTCOMES + ['BU']:
            yrs, D, N, cells = arrays(p, geo, o)
            last = min(yrs.max(), 2023)          # protocol S5: rates of 2021-2023 (2020-2022 where 2023 is missing)
            sel = (yrs > last - 3) & (yrs <= last)
            rate = D[sel].sum(0) / N[sel].sum(0)
            q = pp[(pp.geo == geo) & (pp.variant == 'BSL')].sort_values(['year', 'sex', 'age'])
            tyears = np.array(sorted(q.year.unique()))
            Nt = q.py.to_numpy().reshape(len(tyears), -1)
            for i, t in enumerate(tyears):
                rows.append({'geo': geo, 'outcome': o, 'year': t, 'deaths': float((rate * Nt[i]).sum()),
                             'rate_years': f'{last - 2}-{last}'})
    return pd.DataFrame(rows)


def uk_wpp_sensitivity(p, pp):
    """S11: UK rates 2015-2019 and projection with UN WPP denominators (medium variant)."""
    w = pd.read_csv(RAW / 'unwpp' / 'wpp2024_srb_gbr_age5.csv.gz')
    w = w[(w.ISO3_code == 'GBR') & w.Time.between(2015, 2040)]
    long = w.melt(id_vars=['Time', 'AgeGrpStart'], value_vars=['PopMale', 'PopFemale'], var_name='s', value_name='n')
    long['sex'] = long.s.map({'PopMale': 1, 'PopFemale': 2})
    long['age'] = np.minimum(long.AgeGrpStart, 85)
    wp = long.groupby(['Time', 'sex', 'age']).n.sum() * 1000
    rows = []
    for o in OUTCOMES + ['BU']:
        yrs, D, N, cells = arrays(p, 'UK', o)
        sel = (yrs >= 2015) & (yrs <= 2019)
        Nw = np.vstack([wp.loc[y].reindex(pd.MultiIndex.from_frame(cells)).to_numpy() for y in yrs[sel]])
        rate = D[sel].sum(0) / Nw.sum(0)
        for t in range(2020, 2041):
            nt = wp.loc[t].reindex(pd.MultiIndex.from_frame(cells)).to_numpy()
            rows.append({'geo': 'UK', 'outcome': o, 'year': t, 'deaths_wpp': float((rate * nt).sum())})
    return pd.DataFrame(rows)


def set_w_projection(pp):
    """Malignant outcome, 32 countries: M1 with 2015-2019 rates (Eurostat; UK WHO), BSL and variants."""
    from config import NAMES
    e = pd.read_csv(PROCESSED / 'panel_eurostat.csv.gz')
    e = e[e.sex.isin(['M', 'F'])].assign(sex=lambda x: x.sex.map({'M': 1, 'F': 2}))
    p = load_panel()
    uk = p[p.geo == 'UK'][['geo', 'year', 'sex', 'age', 'MAL', 'pop']]
    e = e[['geo', 'year', 'sex', 'age', 'MAL', 'pop_es']].rename(columns={'pop_es': 'pop'})
    grids = {r.geo: [int(x) for x in r.grid.split()] for r in pd.read_csv(PROCESSED / 'audit' / 'age_grids.csv').itertuples()}
    rows = []
    for g in NAMES:
        gw = [0] + [a for a in grids[g] if a >= 25]
        src = uk if g == 'UK' else e[e.geo == g]
        s = collapse(src[src.year.between(2015, 2019)], gw, value_cols=['MAL', 'pop'])
        assert s.year.nunique() == 5, g
        r = s.groupby(['sex', 'age'])[['MAL', 'pop']].sum()
        rate = r.MAL / r['pop']
        q = collapse(pp[pp.geo == g].rename(columns={'py': 'pop'}), gw, by=('geo', 'variant', 'year', 'sex'), value_cols=['pop'])
        q = q.merge(rate.rename('rate').reset_index(), on=['sex', 'age'], how='left')
        assert q.rate.notna().all(), g
        q['deaths'] = q.rate * q['pop']
        rows.append(q.groupby(['geo', 'variant', 'year'], as_index=False).deaths.sum())
    return pd.concat(rows, ignore_index=True)


def evaluate():
    p, cs, L, bt = setup()
    pp = pd.read_csv(PROCESSED / 'pop_projections.csv.gz')
    bp = MOD / 'bapc_results.csv.gz'
    bapc = pd.read_csv(bp) if bp.exists() else None
    if bapc is None:
        print('WARNING: BAPC results not found; M5 omitted')
    b = backtest(p, bt, bapc)
    b.to_csv(TAB / 'backtest_detail.csv', index=False)
    m, sel = select(b)
    m.to_csv(TAB / 'backtest_summary.csv', index=False)
    sel.to_csv(TAB / 'model_selection.csv', index=False)
    europe_backtest(b).to_csv(TAB / 'backtest_europe.csv', index=False)
    print(m[~m.covid_target].round(2).to_string(index=False))
    print(sel.to_string(index=False))
    pr, agg = project(p, L, bapc, pp)
    pr.to_csv(TAB / 'projections_country.csv', index=False)
    agg.to_csv(TAB / 'projections_aggregate.csv', index=False)
    latest_rates_projection(p, L, pp).to_csv(TAB / 'projections_S5_latest_rates.csv', index=False)
    uk_wpp_sensitivity(p, pp).to_csv(TAB / 'projections_S11_uk_wpp.csv', index=False)
    set_w_projection(pp).to_csv(TAB / 'projections_setW_malignant.csv', index=False)
    a = agg[(agg.geo == 'Set L') & (agg.year.isin([2022, 2030, 2035, 2040])) & (agg.variant == 'BSL')]
    print(a[['outcome', 'model', 'year', 'deaths', 'pi_lo', 'pi_hi']].round(0).to_string(index=False))


if __name__ == '__main__':
    {'prepare': prepare, 'evaluate': evaluate}[sys.argv[1]]()
