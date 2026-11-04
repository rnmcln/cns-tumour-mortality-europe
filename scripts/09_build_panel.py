"""Step 9: analysis panel (protocol v1.0, sections 3-4).

Inputs : data/interim/who_detail.csv.gz (step 8); Eurostat demo_pjan and hlth_cd_aro; WHO population.
Outputs (data/processed):
  panel_who.csv.gz       WHO deaths by country, year (2000-2024), sex and age group, for the elementary
                         groups in config.GROUPS, with Eurostat person-years (pop_es), WHO population
                         (pop_who) and the primary denominator (pop).
  panel_eurostat.csv.gz  Eurostat deaths (C70-C72, C00-C97, R00-R99, all causes) by country, year, sex
                         (M, F; T kept separately) and age group, with Eurostat person-years.
  audit/*.csv            data audit tables (gap-filling, excluded deaths, age grids, denominators).

Rules (protocol section 4)
- Age groups: lower bounds 0, 5, ..., 85, 90, 95. A country's grid ends at 85+ if any of its WHO
  country-years in 2000-2023 (deaths or, for the UK, population) is reported to 85+ only; otherwise 95+.
  Norway: ages 0-14 combined in all years (restricted-format data from 2017).
- Germany 2021-2022: for each elementary group and sex, the difference between the sum of the
  all-ages totals and the sum of the age-specific cells is allocated to the age groups with a
  suppressed cell in proportion to the group's 2018-2020 deaths at those ages.
- Deaths of unknown sex or age are excluded (counts written to the audit).
- Person-years (Eurostat): mean of the 1 January populations of years t and t+1, single years of age
  grouped; population of unknown age distributed in proportion to the known ages.
- UK denominators: ONS mid-year populations as reported to WHO; 2021 and 2022 (missing) and 2018
  (incomplete in the WHO file, about 1.9 million persons short; deviations table, entry 1) interpolated
  geometrically by age and sex between the neighbouring years with data.
"""
import zipfile

import numpy as np
import pandas as pd

from config import INTERIM, PROCESSED, RAW, COUNTRIES, WINDOW, WHO_CODES, GROUPS

AUD = PROCESSED / 'audit'
AUD.mkdir(parents=True, exist_ok=True)
FINE = list(range(0, 100, 5))                    # 0, 5, ..., 95
COL_AGE = {**{f'Deaths{i}': 0 for i in range(2, 7)},
           **{f'Deaths{7 + i}': 5 * (i + 1) for i in range(16)},   # Deaths7 = 5-9 ... Deaths22 = 80-84
           'Deaths23': 85, 'Deaths24': 90, 'Deaths25': 95}


def elementary_group(cause):
    c3 = cause[:3]
    out = []
    if cause == 'AAA':
        return ['ALL']
    if c3 in ('C70', 'C71', 'C72'):
        out.append('MAL')
        if cause in ('C709', 'C719', 'C729'):
            out.append('MAL_UNSPEC')
        if len(cause) == 3:
            out.append('MAL_3CHAR')
    if c3 in ('D32', 'D33', 'D42', 'D43'):
        out.append(c3)
    # subcodes for the analysis of the non-malignant rate component (additional analysis)
    if c3 == 'D43':
        out.append({'D432': 'D43_2', 'D430': 'D43_01', 'D431': 'D43_01'}.get(cause, 'D43_OTH'))
    if cause == 'C719':
        out.append('C71_9')
    if cause[:4] in ('C793',):
        out.append('C79_3')
    if c3 == 'C80':
        out.append('C80')
    if cause == 'R900':
        out.append('R90_0')
    if c3[0] == 'D' and int(c3[1:3]) <= 48:
        out.append('DALL')
    if cause in ('C751', 'C752', 'C753', 'D352', 'D353', 'D354', 'D443', 'D444', 'D445'):
        out.append('SELLAR')
    if c3[0] == 'C' and int(c3[1:3]) <= 97:
        out.append('CALL')
    if c3[0] == 'R':
        out.append('R')
    return out


def who_deaths():
    d = pd.read_csv(INTERIM / 'who_detail.csv.gz', dtype={'List': str, 'Cause': str, 'Frmat': str})
    d = d[d.Year.between(WINDOW[0], 2024)]
    rows = []
    for cause, grps in {c: elementary_group(c) for c in d.Cause.unique()}.items():
        for g in grps:
            rows.append((cause, g))
    m = pd.DataFrame(rows, columns=['Cause', 'grp'])
    d = d.merge(m, on='Cause')
    excl = d[d.Sex == 9].groupby('grp').Deaths1.sum().rename('deaths_unknown_sex')
    unk_age = d[d.Sex.isin([1, 2])].groupby('grp').Deaths26.sum().rename('deaths_unknown_age')
    pd.concat([excl, unk_age], axis=1).fillna(0).to_csv(AUD / 'excluded_unknown_sex_age.csv')
    d = d[d.Sex.isin([1, 2])]
    age_cols = list(COL_AGE)
    long = d.melt(id_vars=['geo', 'Year', 'Sex', 'Frmat', 'List', 'Cause', 'grp', 'Deaths1'],
                  value_vars=age_cols, var_name='col', value_name='deaths')
    long['age'] = long.col.map(COL_AGE)
    long['suppressed'] = long.deaths.isna() & long.Deaths1.notna() & (long.geo == 'DE') & (long.Year >= 2021) \
        & ~((long.Frmat == '01') & long.col.isin(['Deaths24', 'Deaths25']))
    agg = long.groupby(['geo', 'Year', 'Sex', 'grp', 'age'], as_index=False).agg(
        deaths=('deaths', 'sum'), suppressed=('suppressed', 'any'))
    tot = d.groupby(['geo', 'Year', 'Sex', 'grp'], as_index=False).Deaths1.sum()
    fmt = d.groupby(['geo', 'Year']).agg(frmat=('Frmat', 'first'), icd_list=('List', 'first')).reset_index()
    return agg, tot, fmt


def fill_germany(agg, tot):
    """Allocate suppressed deaths, Germany 2021-2022 (protocol section 4, gap-filling)."""
    log = []
    ref = agg[(agg.geo == 'DE') & agg.Year.between(2018, 2020)].groupby(['Sex', 'grp', 'age']).deaths.sum()
    for (y, s, g), sub in agg[(agg.geo == 'DE') & (agg.Year >= 2021)].groupby(['Year', 'Sex', 'grp']):
        t = tot[(tot.geo == 'DE') & (tot.Year == y) & (tot.Sex == s) & (tot.grp == g)].Deaths1.sum()
        diff = t - sub.deaths.sum()
        cand = sub[sub.suppressed]
        if diff < -1e-9:
            raise ValueError(f'negative remainder DE {y} {s} {g}: {diff}')
        if diff == 0:
            continue
        if cand.empty:
            raise ValueError(f'remainder without suppressed cells DE {y} {s} {g}')
        w = np.array([ref.get((s, g, a), 0.0) for a in cand.age], float)
        rule = 'shares 2018-2020'
        if w.sum() == 0:
            w, rule = np.ones(len(cand)), 'equal shares (no reference deaths)'
        alloc = diff * w / w.sum()
        agg.loc[cand.index, 'deaths'] += alloc
        log.append({'year': y, 'sex': s, 'group': g, 'all_ages_total': t, 'sum_of_cells': t - diff,
                    'allocated': diff, 'n_cells': len(cand), 'rule': rule})
    pd.DataFrame(log).to_csv(AUD / 'germany_2021_2022_allocation.csv', index=False)
    # deaths in rows whose all-ages total is itself suppressed (fewer than 3 deaths): count rows only
    return agg


def suppressed_totals():
    d = pd.read_csv(INTERIM / 'who_detail.csv.gz', dtype={'List': str, 'Cause': str, 'Frmat': str})
    d = d[(d.geo == 'DE') & (d.Year >= 2021) & d.Deaths1.isna()]
    d['grp'] = d.Cause.map(lambda c: ','.join(elementary_group(c)))
    out = d.groupby(['Year', 'grp']).size().rename('rows_with_suppressed_total').reset_index()
    out['deaths_lower'] = out.rows_with_suppressed_total * 1
    out['deaths_upper'] = out.rows_with_suppressed_total * 2
    out.to_csv(AUD / 'germany_suppressed_totals.csv', index=False)


def age_grids(fmt, opn):
    """Top age group: 85+ if deaths are reported to 85+ only in any year of the window; otherwise the
    lowest open age of the Eurostat population (rounded down to a multiple of 5) in the window, at most 95."""
    w = fmt[fmt.Year.between(*WINDOW)]
    o = opn[opn.year.between(WINDOW[0], WINDOW[1] + 1)].groupby('geo').open_age.min()
    top = {}
    for g, s in w.groupby('geo'):
        t = 85 if s.frmat.isin(['01', '02']).any() else 95
        if g in o.index:
            t = min(t, int(o[g] // 5 * 5))
        top[g] = t
    top['UK'] = 85                                   # WHO-reported ONS population ends at 85+
    grids = {g: ([0, 15] if g == 'NO' else [0, 5, 10, 15]) + list(range(20, top[g] + 1, 5)) for g in top}
    pd.DataFrame([{'geo': g, 'top_age': top[g], 'grid': ' '.join(map(str, grids[g]))} for g in sorted(grids)]) \
        .to_csv(AUD / 'age_grids.csv', index=False)
    return grids


def to_grid(df, grids, value_cols, keys):
    out = []
    for g, sub in df.groupby('geo'):
        if g not in grids:
            continue
        gr = np.array(grids[g])
        sub = sub.copy()
        sub['age'] = gr[np.searchsorted(gr, sub.age.to_numpy(), side='right') - 1]
        out.append(sub.groupby(keys + ['age'], as_index=False)[value_cols].sum())
    return pd.concat(out, ignore_index=True)


def eurostat_person_years():
    p = pd.read_csv(RAW / 'eurostat' / 'demo_pjan.csv.gz').rename(columns={'geo\\TIME_PERIOD': 'geo'})
    p = p[p.sex.isin(['M', 'F']) & p.geo.isin(COUNTRIES)]
    years = [c for c in p.columns if c.isdigit() and int(c) >= WINDOW[0]]

    def single(a):
        if a == 'Y_LT1':
            return 0
        if a == 'Y_OPEN':
            return -1                               # open age interval, resolved below
        if a.startswith('Y') and a[1:].isdigit():
            return int(a[1:])
        return np.nan
    p['a1'] = p.age.map(single)
    known = p[p.a1.notna()].melt(id_vars=['geo', 'sex', 'a1'], value_vars=years, var_name='year', value_name='n')
    unk = p[p.age == 'UNK'].melt(id_vars=['geo', 'sex'], value_vars=years, var_name='year', value_name='unk')
    known = known.dropna(subset=['n'])
    # the open interval starts one year above the highest single age reported for that country-sex-year
    top = known[known.a1 >= 0].groupby(['geo', 'sex', 'year']).a1.max().rename('open_age') + 1
    known = known.merge(top.reset_index(), on=['geo', 'sex', 'year'])
    known.loc[known.a1 < 0, 'a1'] = known.open_age
    opn = top.reset_index()
    opn['year'] = opn.year.astype(int)
    opn.groupby(['geo', 'year']).open_age.min().reset_index().to_csv(AUD / 'eurostat_population_open_age.csv', index=False)
    known = known.merge(unk, on=['geo', 'sex', 'year'], how='left').fillna({'unk': 0})
    known['n'] = known.n * (1 + known.unk / known.groupby(['geo', 'sex', 'year']).n.transform('sum'))
    known['year'] = known.year.astype(int)
    known['age'] = (known.a1 // 5 * 5).clip(upper=95).astype(int)
    jan = known.groupby(['geo', 'sex', 'year', 'age'], as_index=False).n.sum()
    # the top (open) group of each 1 January population; when the open age differs between 1 January of
    # t and t+1, both are collapsed to the lower of the two top groups before averaging
    jtop = jan.groupby(['geo', 'sex', 'year']).age.max().rename('top').reset_index()
    pair = jtop.merge(jtop.assign(year=jtop.year - 1).rename(columns={'top': 'top_next'}), on=['geo', 'sex', 'year'])
    pair['top_pair'] = pair[['top', 'top_next']].min(axis=1)
    cur = jan.merge(pair[['geo', 'sex', 'year', 'top_pair']], on=['geo', 'sex', 'year'])
    nxt = jan.assign(year=jan.year - 1).merge(pair[['geo', 'sex', 'year', 'top_pair']], on=['geo', 'sex', 'year'])
    for f in (cur, nxt):
        f['age'] = np.minimum(f.age, f.top_pair)
    cur = cur.groupby(['geo', 'sex', 'year', 'age'], as_index=False).n.sum()
    nxt = nxt.groupby(['geo', 'sex', 'year', 'age'], as_index=False).n.sum().rename(columns={'n': 'n_next'})
    py = cur.merge(nxt, on=['geo', 'sex', 'year', 'age'], how='outer', indicator=True)
    assert (py._merge == 'both').all()
    py['pop_es'] = (py.n + py.n_next) / 2
    py['Sex'] = py.sex.map({'M': 1, 'F': 2})
    return py[['geo', 'year', 'Sex', 'age', 'pop_es']], opn


def who_population():
    with zipfile.ZipFile(RAW / 'who_mdb' / 'mort_pop.zip') as z:
        p = pd.read_csv(z.open(z.namelist()[0]), dtype=str)
    p['Country'] = p.Country.astype(int)
    p = p[p.Country.isin(WHO_CODES) & p.Admin1.isna() & p.SubDiv.isna()].copy()
    p['geo'] = p.Country.map(WHO_CODES)
    p['year'] = p.Year.astype(int)
    p = p[p.year.between(WINDOW[0], 2024) & p.Sex.isin(['1', '2'])]
    p['Sex'] = p.Sex.astype(int)
    cols = {f'Pop{i}': a for i, a in zip(range(2, 26), [0] * 5 + list(range(5, 85, 5)) + [85, 90, 95])}
    for c in cols:
        p[c] = pd.to_numeric(p[c], errors='coerce')
    long = p.melt(id_vars=['geo', 'year', 'Sex', 'Frmat'], value_vars=list(cols), var_name='col', value_name='pop_who')
    long['age'] = long.col.map(cols)
    long = long.groupby(['geo', 'year', 'Sex', 'age'], as_index=False).pop_who.sum(min_count=1)
    return long.dropna(subset=['pop_who'])


def interpolate_uk(pw, grid):
    uk = pw[pw.geo == 'UK'].copy()
    gr = np.array(grid)
    uk['age'] = gr[np.searchsorted(gr, uk.age.to_numpy(), side='right') - 1]
    uk = uk.groupby(['geo', 'year', 'Sex', 'age'], as_index=False).pop_who.sum()
    piv = uk.pivot_table(index=['Sex', 'age'], columns='year', values='pop_who')
    log = []
    for y, (y0, y1) in {2018: (2017, 2019), 2021: (2020, 2023), 2022: (2020, 2023)}.items():
        f = (y - y0) / (y1 - y0)
        old = piv[y].sum() if y in piv.columns else np.nan
        piv[y] = piv[y0] * (piv[y1] / piv[y0]) ** f
        log.append({'year': y, 'between': f'{y0}-{y1}', 'who_reported_total': old, 'interpolated_total': piv[y].sum()})
    pd.DataFrame(log).to_csv(AUD / 'uk_population_interpolation.csv', index=False)
    out = piv.stack().rename('pop_uk').reset_index()
    out['geo'] = 'UK'
    return out


def eurostat_deaths(py):
    e = pd.read_csv(RAW / 'eurostat' / 'hlth_cd_aro.csv.gz').rename(columns={'geo\\TIME_PERIOD': 'geo'})
    e = e[e.geo.isin(COUNTRIES) & e.sex.isin(['M', 'F', 'T']) & (e.resid == 'TOT_IN') & (e.unit == 'NR')]
    # ages under 25 are published incompletely in some country-years (single groups blank while the
    # aggregates Y_LT15 / Y_LT25 are given), so ages 0-24 form one group taken from Y_LT25
    amap = {'Y_LT25': 0, 'Y_GE95': 95, **{f'Y{a}-{a + 4}': a for a in range(25, 95, 5)}}
    years = [c for c in e.columns if c.isdigit()]
    tot = e[e.age == 'TOTAL']
    e = e[e.age.isin(amap)].copy()
    e['age_lo'] = e.age.map(amap)
    long = e.melt(id_vars=['geo', 'sex', 'icd10', 'age_lo'], value_vars=years, var_name='year', value_name='deaths')
    long = long.dropna(subset=['deaths'])
    long['year'] = long.year.astype(int)
    out = long.groupby(['geo', 'year', 'sex', 'icd10', 'age_lo'], as_index=False).deaths.sum()
    out = out.pivot_table(index=['geo', 'year', 'sex', 'age_lo'], columns='icd10', values='deaths').reset_index()
    out = out.rename(columns={'age_lo': 'age', 'C70-C72': 'MAL', 'C': 'CALL', 'R': 'R', 'TOTAL': 'ALL',
                              'D00-D48': 'D00_D48'})
    # check: age groups sum to the published total
    tl = tot.melt(id_vars=['geo', 'sex', 'icd10'], value_vars=years, var_name='year', value_name='t').dropna()
    tl['year'] = tl.year.astype(int)
    chk = long.groupby(['geo', 'year', 'sex', 'icd10']).deaths.sum().reset_index().merge(tl, on=['geo', 'year', 'sex', 'icd10'])
    chk['diff'] = chk.t - chk.deaths
    chk[chk['diff'].abs() > 0.5].to_csv(AUD / 'eurostat_age_sum_mismatch.csv', index=False)
    py = py.assign(age=np.where(py.age < 25, 0, py.age)).groupby(['geo', 'year', 'Sex', 'age'], as_index=False).pop_es.sum(min_count=1)
    pye = py.assign(sex=py.Sex.map({1: 'M', 2: 'F'}))[['geo', 'year', 'sex', 'age', 'pop_es']]
    pyt = pye.groupby(['geo', 'year', 'age'], as_index=False).pop_es.sum().assign(sex='T')
    out = out.merge(pd.concat([pye, pyt]), on=['geo', 'year', 'sex', 'age'], how='left')
    return out


def main():
    agg, tot, fmt = who_deaths()
    fmt.to_csv(AUD / 'who_formats_lists.csv', index=False)
    agg = fill_germany(agg, tot)
    suppressed_totals()
    py, opn = eurostat_person_years()
    grids = age_grids(fmt, opn)
    wide = agg.pivot_table(index=['geo', 'Year', 'Sex', 'age'], columns='grp', values='deaths', fill_value=0).reset_index()
    for g in GROUPS:
        if g not in wide.columns:
            wide[g] = 0.0
    wide = wide.rename(columns={'Year': 'year'})
    wide = to_grid(wide, grids, GROUPS, ['geo', 'year', 'Sex'])
    pyg = to_grid(py, grids, ['pop_es'], ['geo', 'year', 'Sex'])
    pw = who_population()
    pwg = to_grid(pw, grids, ['pop_who'], ['geo', 'year', 'Sex'])
    uk = interpolate_uk(pw, grids['UK'])
    panel = wide.merge(pyg, on=['geo', 'year', 'Sex', 'age'], how='left') \
                .merge(pwg, on=['geo', 'year', 'Sex', 'age'], how='left') \
                .merge(uk, on=['geo', 'year', 'Sex', 'age'], how='left')
    panel['pop'] = np.where(panel.geo == 'UK', panel.pop_uk, panel.pop_es)
    panel = panel.merge(fmt.rename(columns={'Year': 'year'}), on=['geo', 'year'], how='left')
    panel = panel.rename(columns={'Sex': 'sex'}).sort_values(['geo', 'year', 'sex', 'age'])
    panel = panel[['geo', 'year', 'sex', 'age', 'frmat', 'icd_list'] + GROUPS + ['pop', 'pop_es', 'pop_who', 'pop_uk']]
    panel.to_csv(PROCESSED / 'panel_who.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
    # denominator availability by country-year
    av = panel.groupby(['geo', 'year']).agg(deaths_all=('ALL', 'sum'), pop=('pop', lambda s: s.notna().all()))
    av.reset_index().to_csv(AUD / 'panel_availability.csv', index=False)
    eu = eurostat_deaths(py)
    eu.to_csv(PROCESSED / 'panel_eurostat.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
    print('panel rows', len(panel), '| countries', panel.geo.nunique(), '| years', panel.year.min(), panel.year.max())
    print('missing primary pop rows:', panel['pop'].isna().sum())


if __name__ == '__main__':
    main()
