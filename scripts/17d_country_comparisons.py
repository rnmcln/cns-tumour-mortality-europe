"""Step 17d: country comparisons on malignant-only and all-CNS definitions (deviations table, entry 20).

For the Set L countries with country-level results, age-standardised death rates (2013 ESP, both sexes, grid
0-14, 15-19, ..., 85+) are calculated for 2003-2005 and 2020-2022 from pooled deaths and person-years for
the malignant outcome (C70-C72), the all-CNS outcome (C70-C72, D32-D33, D42-D43) and the narrow all-CNS
outcome (C70-C72, D33, deviations table, entry 17). Countries are ranked from highest to lowest rate on each definition.
Uncertainty: 2,000 Poisson draws of the death counts in every country, period and age group (seed from
config), from which 95% intervals of each country's rank and of its rank difference are taken.
Direction of change: percentage change in each age-standardised rate between the two periods, with 95%
intervals from the same draws.
Outputs: output/tables/comparison_country_ranks.csv, comparison_rank_summary.csv,
         comparison_trend_direction.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from config import OUTPUT, ROOT, SEED, BASE_L, END_L
from lib.data import load_panel, collapse, esp_weights

TAB = OUTPUT / 'tables'
GRID = [0, 15] + list(range(20, 90, 5))
DEF = {'MAL': ['MAL'], 'COMP': ['MAL', 'D32', 'D33', 'D42', 'D43'], 'NARROW': ['MAL', 'D33', 'D43']}
NSIM = 2000


def main():
    p = load_panel()
    p = p[p['pop'].notna()]
    cs = pd.read_csv(ROOT / 'protocol' / 'country_sets.csv')
    geos = cs[cs.set_L & cs.country_level_reporting].geo.tolist()
    w = esp_weights(GRID)
    w = (w / w.sum()).to_numpy()
    elem = ['MAL', 'D32', 'D33', 'D42', 'D43']
    rng = np.random.default_rng(SEED + 40)
    rates, sims = {}, {}
    for per, (a, b) in [('0', BASE_L), ('1', END_L)]:
        q = collapse(p[p.geo.isin(geos) & p.year.between(a, b)], GRID, by=('geo',), value_cols=elem + ['pop'])
        for g in geos:
            s = q[q.geo == g].set_index('age').reindex(sorted(GRID))
            n = s['pop'].to_numpy()
            draws = {e: rng.poisson(np.broadcast_to(s[e].to_numpy(), (NSIM, len(n)))) for e in elem}
            for d, parts in DEF.items():
                dd = sum(s[e].to_numpy() for e in parts)
                rates[(d, per, g)] = float((dd / n * w).sum() * 1e5)
                sims[(d, per, g)] = (sum(draws[e] for e in parts) / n * w).sum(1) * 1e5
    rows, summ = [], []
    for per in ['0', '1']:
        rk = {}
        for d in DEF:
            v = pd.Series({g: rates[(d, per, g)] for g in geos})
            rk[d] = v.rank(ascending=False, method='min')
            simv = np.column_stack([sims[(d, per, g)] for g in geos])
            rk[d + '_sim'] = (-simv).argsort(1).argsort(1) + 1          # rank 1 = highest rate
        for g_i, g in enumerate(geos):
            diff_sim = rk['MAL_sim'][:, g_i] - rk['COMP_sim'][:, g_i]
            rows.append({'geo': g, 'period': per, 'asr_MAL': rates[('MAL', per, g)], 'asr_COMP': rates[('COMP', per, g)],
                         'asr_NARROW': rates[('NARROW', per, g)], 'rank_MAL': int(rk['MAL'][g]), 'rank_COMP': int(rk['COMP'][g]),
                         'rank_NARROW': int(rk['NARROW'][g]), 'rank_diff': int(rk['MAL'][g] - rk['COMP'][g]),
                         'rank_diff_lo': float(np.percentile(diff_sim, 2.5)), 'rank_diff_hi': float(np.percentile(diff_sim, 97.5))})
        t = pd.DataFrame([r for r in rows if r['period'] == per])
        ci_excl0 = ((t.rank_diff_lo > 0) | (t.rank_diff_hi < 0))
        summ.append({'period': per, 'n_countries': len(t),
                     'spearman_MAL_COMP': spearmanr(t.asr_MAL, t.asr_COMP)[0],
                     'spearman_MAL_NARROW': spearmanr(t.asr_MAL, t.asr_NARROW)[0],
                     'mean_abs_rank_diff': t.rank_diff.abs().mean(), 'max_abs_rank_diff': int(t.rank_diff.abs().max()),
                     'n_shift_ge5': int((t.rank_diff.abs() >= 5).sum()), 'n_shift_ci_excl0': int(ci_excl0.sum()),
                     'n_shift_ge5_ci_excl0': int(((t.rank_diff.abs() >= 5) & ci_excl0).sum())})
    pd.DataFrame(rows).to_csv(TAB / 'comparison_country_ranks.csv', index=False)
    pd.DataFrame(summ).to_csv(TAB / 'comparison_rank_summary.csv', index=False)
    tr = []
    for g in geos:
        r = {'geo': g}
        for d in ['MAL', 'COMP']:
            ch = 100 * (sims[(d, '1', g)] / sims[(d, '0', g)] - 1)
            r[f'pct_{d}'] = 100 * (rates[(d, '1', g)] / rates[(d, '0', g)] - 1)
            r[f'pct_{d}_lo'], r[f'pct_{d}_hi'] = np.percentile(ch, [2.5, 97.5])
        r['opposite'] = bool(np.sign(r['pct_MAL']) != np.sign(r['pct_COMP']))
        r['opposite_both_excl0'] = bool(r['opposite'] and (r['pct_MAL_lo'] > 0 or r['pct_MAL_hi'] < 0)
                                        and (r['pct_COMP_lo'] > 0 or r['pct_COMP_hi'] < 0))
        tr.append(r)
    tr = pd.DataFrame(tr)
    tr.to_csv(TAB / 'comparison_trend_direction.csv', index=False)
    print(pd.DataFrame(summ).round(2).to_string(index=False))
    print(tr.round(1).sort_values('pct_MAL').to_string(index=False))


if __name__ == '__main__':
    main()
