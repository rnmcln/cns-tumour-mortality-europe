"""Step 10: coding changes (protocol v1.0, section 5.2).

Documented changes (2011 onwards): data/interim/coding_change_points.csv (step 6, protocol appendix A).
Data-identified candidates before 2011, in countries with at least 100 C70-C72 deaths per year
(country_level_reporting in protocol/country_sets.csv) and WHO data for both outcomes:
  year t (t <= 2010) is a candidate if
  (a) |mean_{t..t+2} log(BU/MAL) - mean_{t-3..t-1} log(BU/MAL)| > 0.3, where BU/MAL is the annual ratio of
      benign and uncertain-behaviour (D32-D33, D42-D43) to malignant (C70-C72) CNS deaths, and
  (b) the age-standardised death rates of the malignant and of the benign or uncertain series change in
      opposite directions between the two windows (rates from pooled deaths and person-years per window).
  All six years must be observed. Operational clarification (not stated in the protocol): consecutive
  candidate years are one change, dated at the year with the largest difference in (a).
Outputs: output/tables/coding_rule_statistics.csv (every evaluated country-year),
         data/processed/coding_changes.csv (documented changes and data-identified candidates).
"""
import numpy as np
import pandas as pd

from config import INTERIM, PROCESSED, OUTPUT, ROOT, RAW
from lib.data import load_panel, collapse, asr

TAB = OUTPUT / 'tables'
TAB.mkdir(parents=True, exist_ok=True)
GRID = [0, 15] + list(range(20, 90, 5))           # common grid 0-14, 15-19, ..., 85+


def window_asr(sub, years, col):
    w = sub[sub.year.isin(years)].groupby('age', as_index=False)[[col, 'pop']].sum()
    return asr(w, col, 'pop', GRID)


def main():
    p = load_panel()
    p = p[p['pop'].notna()]
    cs = pd.read_csv(ROOT / 'protocol' / 'country_sets.csv')
    eligible = cs[cs.country_level_reporting].geo.tolist()
    q = collapse(p, GRID, by=('geo', 'year'), value_cols=['MAL', 'BU', 'pop'])
    ann = q.groupby(['geo', 'year'], as_index=False)[['MAL', 'BU']].sum()
    ann['lr'] = np.log(ann.BU / ann.MAL)
    rows = []
    for g in eligible:
        a = ann[ann.geo == g].set_index('year')
        sub = q[q.geo == g]
        for t in range(2003, 2011):
            pre, post = list(range(t - 3, t)), list(range(t, t + 3))
            if not all(y in a.index for y in pre + post):
                continue
            if not np.isfinite(a.loc[pre + post, 'lr']).all():
                continue
            d_lr = a.loc[post, 'lr'].mean() - a.loc[pre, 'lr'].mean()
            m0, m1 = window_asr(sub, pre, 'MAL'), window_asr(sub, post, 'MAL')
            b0, b1 = window_asr(sub, pre, 'BU'), window_asr(sub, post, 'BU')
            opp = np.sign(m1 - m0) * np.sign(b1 - b0) < 0
            rows.append({'geo': g, 'year': t, 'diff_mean_log_ratio': d_lr, 'asr_mal_pre': m0, 'asr_mal_post': m1,
                         'asr_bu_pre': b0, 'asr_bu_post': b1, 'opposite_directions': bool(opp),
                         'candidate': bool(abs(d_lr) > 0.3 and opp)})
    st = pd.DataFrame(rows)
    st.to_csv(TAB / 'coding_rule_statistics.csv', index=False)
    cand = []
    for g, s in st[st.candidate].groupby('geo'):
        s = s.sort_values('year')
        run = (s.year.diff() != 1).cumsum()
        for _, r in s.groupby(run):
            best = r.loc[r.diff_mean_log_ratio.abs().idxmax()]
            cand.append({'geo': g, 'year': int(best.year), 'source': 'data rule',
                         'detail': f'difference in mean log ratio {best.diff_mean_log_ratio:+.2f}; '
                                   f'candidate years {r.year.min()}-{r.year.max()}'})
    doc = pd.read_csv(INTERIM / 'coding_change_points.csv')
    doc = pd.DataFrame({'geo': doc.geo, 'year': doc.year, 'source': 'documented',
                        'detail': doc['from'] + ' -> ' + doc['to']})
    # changes documented outside the Eurostat metadata (deviations table, entry 13): UK 2014 (ONS); Germany undated
    add = pd.read_csv(RAW / 'metadata' / 'additional_coding_changes.csv')
    dated = add[add.year.notna()]
    doc = pd.concat([doc, pd.DataFrame({'geo': dated.geo, 'year': dated.year.astype(int), 'source': 'documented',
                                        'detail': dated['from'] + ' -> ' + dated['to'] + ' (' + dated['source'] + ')'})])
    add[add.year.isna()][['geo', 'source', 'note']].to_csv(PROCESSED / 'coding_changes_undated.csv', index=False)
    out = pd.concat([doc, pd.DataFrame(cand)], ignore_index=True).sort_values(['geo', 'year'])
    out.to_csv(PROCESSED / 'coding_changes.csv', index=False)
    print(out.to_string(index=False))
    print('evaluated country-years:', len(st), '| countries:', st.geo.nunique())


if __name__ == '__main__':
    main()
