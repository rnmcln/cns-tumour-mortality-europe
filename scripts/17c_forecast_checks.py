"""Step 17c: checks of the expected-death and projection analyses (additional analyses).

1. Constant-rate expectations: for each origin t (2010, 2013, 2016, 2019), age-sex-specific rates of
   t-2..t applied to observed person-years of t+1..t+3 (and 2023 for origin 2019), Set L countries with
   data in all years concerned. Pre-pandemic origins give the empirical distribution of observed to
   expected ratios against which 2020-2022 is compared.
2. Hierarchical-model pseudo-origins: observed to expected ratios from 12_trend_model.R fitted to 2010,
   2013 and 2016 (files trend_obs_exp_val<year>.csv), collected in one table.
3. Aggregate backtest errors of the projection models per origin and horizon (backtest_europe.csv) and an
   empirical calibration of the constant-rate projection: the log ratio of predicted to observed aggregate
   deaths at horizons 5 and 10 (targets outside 2020-2022) is modelled as b*h + e, var(e) = s2*h; b and s2
   are estimated from these points and applied to the 2040 projection (h = 21) and 2030 (h = 11). With
   five points and extrapolation beyond the evaluated horizon, the result is indicative only.
Outputs: output/tables/forecast_checks_*.csv
"""
import numpy as np
import pandas as pd

from config import OUTPUT, ROOT
from lib.data import load_panel

TAB = OUTPUT / 'tables'
BU = ['D32', 'D33', 'D42', 'D43']


def main():
    p = load_panel()
    p = p[p['pop'].notna()].copy()
    p['BU'] = p[BU].sum(axis=1)
    p['COMP'] = p.MAL + p.BU
    cs = pd.read_csv(ROOT / 'protocol' / 'country_sets.csv')
    L = cs[cs.set_L].geo.tolist()
    p = p[p.geo.isin(L)]
    # ---- 1. constant-rate expectations ------------------------------------------------------------
    rows = []
    for t in [2010, 2013, 2016, 2019]:
        base = [t - 2, t - 1, t]
        targets = [t + 1, t + 2, t + 3] + ([2023] if t == 2019 else [])
        for tg in targets:
            for blockname, tl in [('single', [tg])] + ([('3y', [t + 1, t + 2, t + 3])] if tg == t + 1 else []):
                need = base + tl
                ok = [g for g in L if set(need) <= set(p[p.geo == g].year)]
                b = p[p.geo.isin(ok) & p.year.isin(base)].groupby(['geo', 'sex', 'age'])[['MAL', 'BU', 'COMP', 'pop']].sum()
                f = p[p.geo.isin(ok) & p.year.isin(tl)].groupby(['geo', 'sex', 'age'])[['MAL', 'BU', 'COMP', 'pop']].sum()
                for o in ['MAL', 'BU', 'COMP']:
                    exp_ = (b[o] / b['pop'] * f['pop']).sum()
                    rows.append({'origin': t, 'target': tg if blockname == 'single' else f'{t + 1}-{t + 3}', 'outcome': o,
                                 'n_countries': len(ok), 'observed': f[o].sum(), 'expected': exp_, 'ratio': f[o].sum() / exp_})
    cr = pd.DataFrame(rows)
    cr.to_csv(TAB / 'forecast_checks_constant_rate_oe.csv', index=False)
    # ---- 2. hierarchical pseudo-origins ----------------------------------------------------------------
    hv = []
    for fe in [2010, 2013, 2016]:
        f = TAB / f'trend_obs_exp_val{fe}.csv'
        if f.exists():
            x = pd.read_csv(f)
            x = x[x.geo == 'Set L'].assign(fit_end=fe)
            ex = TAB / f'trend_excluded_val{fe}.csv'
            x['excluded'] = ' '.join(pd.read_csv(ex).iloc[:, 0].astype(str)) if ex.exists() and ex.stat().st_size > 30 else ''
            hv.append(x)
    main_oe = pd.read_csv(TAB / 'trend_obs_exp.csv')
    hv.append(main_oe[main_oe.geo == 'Set L'].assign(fit_end=2019, excluded=''))
    hv = pd.concat(hv, ignore_index=True)
    hv.to_csv(TAB / 'forecast_checks_hierarchical_oe.csv', index=False)
    # ---- 3. aggregate backtest and calibration -----------------------------------------------------------
    be = pd.read_csv(TAB / 'backtest_europe.csv')
    be['covid_target'] = be.target.between(2020, 2022)
    be.to_csv(TAB / 'forecast_checks_backtest_aggregate.csv', index=False)
    agg = pd.read_csv(TAB / 'projections_aggregate.csv')
    cal = []
    for o in ['MAL', 'COMP']:
        e = be[(be.outcome == o) & (be.model == 'M1 constant') & ~be.covid_target]
        h = e.horizon.to_numpy(float)
        y = np.log(e.predicted / e.observed).to_numpy()
        b = np.sum(y) / np.sum(h)                          # WLS with var proportional to h: sum(y)/sum(h)
        s2 = np.sum((y - b * h) ** 2 / h) / (len(y) - 1)
        for yr in [2030, 2040]:
            H = yr - 2019
            r = agg[(agg.geo == 'Set L') & (agg.outcome == o) & (agg.model == 'M1 constant') & (agg.variant == 'BSL') & (agg.year == yr)].iloc[0]
            adj = r.deaths * np.exp(-b * H)
            cal.append({'outcome': o, 'year': yr, 'horizon': H, 'n_points': len(y), 'bias_per_year_pct': 100 * (np.exp(b) - 1),
                        'projection': r.deaths, 'adjusted': adj, 'lo': adj * np.exp(-1.96 * np.sqrt(s2 * H)),
                        'hi': adj * np.exp(1.96 * np.sqrt(s2 * H)),
                        'unadj_lo': r.deaths * np.exp(-1.96 * np.sqrt(s2 * H)), 'unadj_hi': r.deaths * np.exp(1.96 * np.sqrt(s2 * H))})
    cal = pd.DataFrame(cal)
    cal.to_csv(TAB / 'forecast_checks_calibration.csv', index=False)
    print(cr[cr.target.astype(str).str.contains('-') | (cr.target.astype(str) == '2023')].round(3).to_string(index=False))
    print(hv[['fit_end', 'outcome', 'year', 'observed', 'expected', 'pi_lo', 'pi_hi', 'ratio', 'excluded']].round(3).to_string(index=False))
    print(be[be.model == 'M1 constant'][['outcome', 'origin', 'horizon', 'target', 'pe', 'pe_85', 'n_countries']].round(2).to_string(index=False))
    print(cal.round(1).to_string(index=False))


if __name__ == '__main__':
    main()
