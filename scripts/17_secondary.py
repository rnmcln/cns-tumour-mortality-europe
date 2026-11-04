"""Step 17: secondary analyses and remaining sensitivity tables (protocol v1.0, sections 3, 5.2, 5.5, 5.6).

Outputs (output/tables):
  bu_by_country_period.csv      mean annual deaths of benign and uncertain-behaviour CNS tumours (D32, D33,
                                D42, deviations table, entry 17) and of C70-C72, ratio BU/MAL, share at 75+ and 85+, Set L, 2003-2005,
                                2011-2013 and 2020-2022
  bu_by_age_sex.csv             Set L totals by broad age group, sex and period
  ratio_dispersion.csv          between-country distribution of the ratio BU/MAL by period
  cns_share_of_cancer.csv       MAL and COMP as a share of C00-C97 deaths by year (Set L)
  decomp_all_cancer.csv         decomposition of C00-C97 deaths, Set L, 2003-2005 to 2020-2022 (reference)
  certification_indicators.csv  share of C70.9/C71.9/C72.9 among C70-C72 deaths (country-years with all
                                C70-C72 deaths at four characters) and share of deaths at 75+ assigned to
                                R00-R99, by broad age group and period
  who_vs_eurostat.csv           S6: C70-C72 deaths, WHO against Eurostat, by country and year (2011-2023)
  step_direction.csv            step terms of the trend model (MAL and BU) and the opposite-direction rule
"""
import numpy as np
import pandas as pd

from config import OUTPUT, PROCESSED, ROOT, REGION, BASE_L, END_L
from lib.data import load_panel
from lib.decomp import decompose

TAB = OUTPUT / 'tables'
PERIODS = {'2003-2005': (2003, 2005), '2011-2013': (2011, 2013), '2020-2022': (2020, 2022)}
BROAD = [(0, 39, '0-39'), (40, 64, '40-64'), (65, 74, '65-74'), (75, 84, '75-84'), (85, 200, '85+')]


def broad(a):
    for lo, hi, lab in BROAD:
        if lo <= a <= hi:
            return lab


def main():
    p = load_panel()
    cs = pd.read_csv(ROOT / 'protocol' / 'country_sets.csv')
    L = cs[cs.set_L].geo.tolist()
    pl = p[p.geo.isin(L)].copy()
    pl['broad'] = pl.age.map(broad)

    rows = []
    for lab, (a, b) in PERIODS.items():
        s = pl[pl.year.between(a, b)]
        ny = s.groupby('geo').year.nunique()
        for g in L:
            t = s[s.geo == g]
            if ny.get(g, 0) != b - a + 1:
                continue
            m = t[['MAL', 'D32', 'D33', 'D42', 'D43', 'BU', 'COMP']].sum() / (b - a + 1)
            rows.append({'geo': g, 'period': lab, **m.round(3).to_dict(), 'ratio_BU_MAL': m.BU / m.MAL,
                         'share75_MAL': t[t.age >= 75].MAL.sum() / t.MAL.sum(), 'share85_MAL': t[t.age >= 85].MAL.sum() / t.MAL.sum(),
                         'share75_BU': t[t.age >= 75].BU.sum() / t.BU.sum(), 'share85_BU': t[t.age >= 85].BU.sum() / t.BU.sum(),
                         'share75_COMP': t[t.age >= 75].COMP.sum() / t.COMP.sum(), 'share85_COMP': t[t.age >= 85].COMP.sum() / t.COMP.sum()})
    bu = pd.DataFrame(rows)
    # age-standardised rates per period (pooled deaths and person-years; grid 0-14, 15-19, ..., 85+; ESP 2013)
    from lib.data import collapse, asr
    grid = [0, 15] + list(range(20, 90, 5))
    qa = collapse(pl, grid, value_cols=['MAL', 'BU', 'COMP', 'pop'])
    for o in ['MAL', 'BU', 'COMP']:
        bu[f'asr_{o}'] = [asr(qa[(qa.geo == r.geo) & qa.year.between(*PERIODS[r.period])].groupby('age', as_index=False)[[o, 'pop']].sum(),
                              o, 'pop', grid) for r in bu.itertuples()]
    eu_rows = []
    for lab, (a, b) in PERIODS.items():
        full = bu[bu.period == lab].geo.tolist()
        t = qa[qa.geo.isin(full) & qa.year.between(a, b)].groupby('age', as_index=False)[['MAL', 'BU', 'COMP', 'pop']].sum()
        eu_rows.append({'period': lab, 'n_countries': len(full), **{f'asr_{o}': asr(t, o, 'pop', grid) for o in ['MAL', 'BU', 'COMP']}})
    pd.DataFrame(eu_rows).to_csv(TAB / 'asr_by_period_setL.csv', index=False)
    bu.to_csv(TAB / 'bu_by_country_period.csv', index=False)
    disp = bu[bu.geo.isin(cs[cs.country_level_reporting].geo)].groupby('period').ratio_BU_MAL.describe(percentiles=[.25, .5, .75])
    disp.to_csv(TAB / 'ratio_dispersion.csv')

    rows = []
    for lab, (a, b) in PERIODS.items():
        s = pl[pl.year.between(a, b)]
        full = [g for g in L if s[s.geo == g].year.nunique() == b - a + 1]
        s = s[s.geo.isin(full)]
        t = s.groupby(['sex', 'broad'])[['MAL', 'BU', 'D32', 'D33', 'D42', 'D43', 'COMP', 'pop']].sum() / (b - a + 1)
        t['rate_BU_per100k'] = t.BU / t['pop'] * 1e5
        t['rate_MAL_per100k'] = t.MAL / t['pop'] * 1e5
        rows.append(t.reset_index().assign(period=lab, n_countries=len(full)))
    pd.concat(rows).to_csv(TAB / 'bu_by_age_sex.csv', index=False)

    yr = pl[pl.year.between(2003, 2022)].groupby('year')[['MAL', 'COMP', 'BU', 'CALL', 'ALL']].sum()
    yr['MAL_share_of_CALL'] = yr.MAL / yr.CALL
    yr['COMP_share_of_CALL_plus_BU'] = yr.COMP / (yr.CALL + yr.BU)
    yr.reset_index().to_csv(TAB / 'cns_share_of_cancer.csv', index=False)

    rows = []
    for g in L:
        def per(y0, y1):
            s = pl[(pl.geo == g) & pl.year.between(y0, y1)].sort_values(['year', 'sex', 'age'])
            ny = y1 - y0 + 1
            nc = len(s) // ny
            return (s.CALL.to_numpy(float).reshape(ny, nc).mean(0), s['pop'].to_numpy(float).reshape(ny, nc).mean(0))
        d0, n0 = per(*BASE_L)
        d1, n1 = per(*END_L)
        r = decompose(d0, n0, d1, n1)
        rows.append({'geo': g, **{k: float(v) for k, v in r.items()}})
    dc = pd.DataFrame(rows)
    dc = pd.concat([dc, dc.drop(columns='geo').sum().to_frame().T.assign(geo='Europe')], ignore_index=True)
    dc.to_csv(TAB / 'decomp_all_cancer.csv', index=False)

    # certification indicators
    fmt = pd.read_csv(PROCESSED / 'audit' / 'who_formats_lists.csv', dtype={'icd_list': str})
    four = pl.groupby(['geo', 'year']).MAL_3CHAR.sum().rename('mal3').reset_index().merge(
        fmt.rename(columns={'Year': 'year'}), on=['geo', 'year'])
    four = four[four.icd_list.isin(['104', '10M']) & (four.mal3 == 0)][['geo', 'year']]
    q = pl.merge(four.assign(four_char=True), on=['geo', 'year'], how='left')
    rows = []
    for lab, (a, b) in PERIODS.items():
        s = q[q.year.between(a, b)]
        sf = s[s.four_char.fillna(False).astype(bool)]
        for br, t in sf.groupby('broad'):
            rows.append({'indicator': 'unspecified subcode share of C70-C72', 'period': lab, 'broad_age': br,
                         'numerator': t.MAL_UNSPEC.sum(), 'denominator': t.MAL.sum(), 'n_countries': t.geo.nunique()})
        rows.append({'indicator': 'unspecified subcode share of C70-C72', 'period': lab, 'broad_age': 'all',
                     'numerator': sf.MAL_UNSPEC.sum(), 'denominator': sf.MAL.sum(), 'n_countries': sf.geo.nunique()})
        full = [g for g in L if s[s.geo == g].year.nunique() == b - a + 1]
        o = s[s.geo.isin(full) & (s.age >= 75)]
        for br, t in o.groupby('broad'):
            rows.append({'indicator': 'R00-R99 share of all deaths', 'period': lab, 'broad_age': br,
                         'numerator': t.R.sum(), 'denominator': t.ALL.sum(), 'n_countries': t.geo.nunique()})
        rows.append({'indicator': 'R00-R99 share of all deaths', 'period': lab, 'broad_age': '75+',
                     'numerator': o.R.sum(), 'denominator': o.ALL.sum(), 'n_countries': o.geo.nunique()})
    # unspecified subcodes in the countries coded at four characters in every year of both 2003-2005 and 2020-2022
    fc = q[q.four_char.fillna(False).astype(bool)]
    com = [g for g in L if all(fc[(fc.geo == g)].year.isin(range(a, b + 1)).sum() > 0 and
                               set(range(a, b + 1)) <= set(fc[fc.geo == g].year) for a, b in [(2003, 2005), (2020, 2022)])]
    for lab, (a, b) in PERIODS.items():
        t = fc[fc.geo.isin(com) & fc.year.between(a, b)]
        rows.append({'indicator': 'unspecified subcode share of C70-C72, common countries', 'period': lab, 'broad_age': 'all',
                     'numerator': t.MAL_UNSPEC.sum(), 'denominator': t.MAL.sum(), 'n_countries': t.geo.nunique()})
    ci = pd.DataFrame(rows)
    ci['share'] = ci.numerator / ci.denominator
    ci.to_csv(TAB / 'certification_indicators.csv', index=False)
    # by country, unspecified share and R share at 75+, per period
    rows = []
    for lab, (a, b) in PERIODS.items():
        s = q[q.year.between(a, b)]
        for g, t in s.groupby('geo'):
            tf = t[t.four_char.fillna(False).astype(bool)]
            rows.append({'geo': g, 'period': lab, 'unspec_share': tf.MAL_UNSPEC.sum() / tf.MAL.sum() if len(tf) else np.nan,
                         'unspec_years': tf.year.nunique(), 'R_share_75plus': t[t.age >= 75].R.sum() / t[t.age >= 75].ALL.sum()})
    pd.DataFrame(rows).to_csv(TAB / 'certification_by_country.csv', index=False)

    # S6 WHO vs Eurostat counts
    e = pd.read_csv(PROCESSED / 'panel_eurostat.csv.gz')
    et = e[e.sex.isin(['M', 'F'])].groupby(['geo', 'year']).MAL.sum().rename('eurostat')
    wt = p.groupby(['geo', 'year']).MAL.sum().rename('who')
    c = pd.concat([wt, et], axis=1, join='inner').reset_index()
    c = c[c.year.between(2011, 2023)]
    c['diff'] = c.who - c.eurostat
    c['rel_diff_pct'] = 100 * c['diff'] / c.eurostat
    c.to_csv(TAB / 'who_vs_eurostat.csv', index=False)

    # contribution of the years of coding change to the chained coding component (year t-1 to t)
    chn = pd.read_csv(TAB / 'decomp_chained.csv')
    cod = chn[(chn.outcome == 'CODING') & chn.geo.isin(L)].pivot(index='year', columns='geo', values='rates')
    chg_all = pd.read_csv(PROCESSED / 'coding_changes.csv')
    rows = []
    for r in chg_all.itertuples():
        if r.geo in L and 2004 <= r.year <= 2022:
            rows.append({'geo': r.geo, 'year': r.year, 'source': r.source,
                         'coding_contribution_in_change_year': cod.loc[r.year, r.geo] - cod.loc[r.year - 1, r.geo]})
    cs_ = pd.DataFrame(rows)
    cs_.to_csv(TAB / 'coding_change_year_contributions.csv', index=False)
    pd.DataFrame({'quantity': ['chained coding component 2003-2022', 'sum over change years (all)', 'sum over change years (documented)',
                               'sum over change years (data rule)'],
                  'deaths': [cod.loc[2022].sum(), cs_.coding_contribution_in_change_year.sum(),
                             cs_[cs_.source == 'documented'].coding_contribution_in_change_year.sum(),
                             cs_[cs_.source == 'data rule'].coding_contribution_in_change_year.sum()]}) \
        .to_csv(TAB / 'coding_change_year_summary.csv', index=False)

    # opposite-direction rule on the step terms
    st = pd.read_csv(TAB / 'trend_steps.csv')
    w = st[st.outcome.isin(['MAL', 'BU'])].pivot_table(index=['geo', 'year'], columns='outcome', values=['rr', 'lo', 'hi'])
    w.columns = [f'{a}_{b}' for a, b in w.columns]
    w = w.reset_index()
    w['opposite_direction'] = np.sign(np.log(w.rr_MAL)) * np.sign(np.log(w.rr_BU)) < 0
    w['both_cis_exclude_1'] = ((w.lo_MAL > 1) | (w.hi_MAL < 1)) & ((w.lo_BU > 1) | (w.hi_BU < 1))
    ch = pd.read_csv(PROCESSED / 'coding_changes.csv')
    w = w.merge(ch, on=['geo', 'year'], how='left')
    w.to_csv(TAB / 'step_direction.csv', index=False)
    print(bu[bu.period != '2011-2013'].groupby('period')[['MAL', 'BU']].sum())
    print(disp.round(3))
    print(ci.round(4).to_string(index=False))
    print(c.groupby('geo').rel_diff_pct.agg(['min', 'max']).round(2).to_string())
    print(w[['geo', 'year', 'rr_MAL', 'rr_BU', 'opposite_direction', 'both_cis_exclude_1', 'source']].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
