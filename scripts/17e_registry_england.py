"""Step 17e: cancer registry incidence by tumour behaviour, England (deviations table, entry 21).

Question: did the decline in deaths recorded as CNS tumours of uncertain or unknown behaviour (deviations table, entry 17) reflect fewer
such tumours being diagnosed, or a change in how deaths are certified? If D43 registrations did not fall while D43
deaths fell, and malignant (C71) deaths rose relative to C71 registrations, the mortality pattern is not explained
by incidence.

Inputs
- ONS Cancer registration statistics, England, annual releases for diagnosis years 2009-2017
  (data/raw/registry_england/ons_<year>.xls[x]): registrations by ICD-10 code (third and fourth character), sex
  and five-year age group (Table 1 or 'Table 1 Data'; fourth-character table: Table 6, or Table 8 in 2015-2016).
- ONS mid-year population estimates, MYEB1 detailed time series 2001-2020 (ons_pop_detailed_2001_2020.zip),
  England (registrations) and England and Wales (deaths).
- WHO Mortality Database, country code 4310 (United Kingdom, England and Wales): deaths by underlying cause
  (fourth character), sex and age, 2001-2020 (same raw files as step 8).
Registry and death data refer to England and to England and Wales respectively (Wales: about 5% of the population).

Outputs (output/tables)
- registry_england_counts.csv: registrations (England) and deaths (England and Wales) by year, sex, age group, code group
- registry_england_asr.csv: age-standardised rates (2013 ESP, both sexes, grid 0-14, 15-19, ..., 85+) by year
- registry_england_trends.csv: average annual percentage change 2009-2017 from Poisson regression of counts on
  calendar year with age-sex terms and log population offset; 95% CI with quasi-Poisson (Pearson) scale
- registry_england_ratios.csv: mortality to incidence ratios of age-standardised rates and subcode shares by year
"""
import io
import re
import zipfile

import numpy as np
import pandas as pd
import statsmodels.api as sm

from config import RAW, OUTPUT
from lib.data import esp_weights

REG = RAW / 'registry_england'
TAB = OUTPUT / 'tables'
GRID = [0, 15] + list(range(20, 90, 5))
YEARS = list(range(2009, 2018))
CODES3 = ['C70', 'C71', 'C72', 'D32', 'D33', 'D42', 'D43']
CODES4 = ['C71.9', 'D43.0', 'D43.1', 'D43.2']
GROUPS = {'MAL': ['C70', 'C71', 'C72'], 'C71': ['C71'], 'C71.9': ['C71.9'], 'D43': ['D43'], 'D43.2': ['D43.2'],
          'D43.0-1': ['D43.0', 'D43.1'], 'D32': ['D32'], 'D33': ['D33'], 'D42': ['D42'],
          'BU': ['D32', 'D33', 'D42', 'D43'], 'COMP': ['C70', 'C71', 'C72', 'D32', 'D33', 'D42', 'D43']}
# lower bounds of the age groups in the ONS tables, in column order after 'All ages'
AGES19 = [0, 1] + list(range(5, 90, 5))                 # Under 1, 1-4, 5-9, ..., 80-84, 85+
AGES20 = [0, 1] + list(range(5, 95, 5))                 # ..., 85-89, 90+


def to_grid(a):
    g = np.asarray(GRID)
    return int(g[np.searchsorted(g, a, side='right') - 1])


def num(v):
    if isinstance(v, str):
        v = v.strip()
        if v in ('-', '', 'u', ':', '*'):
            return 0.0
        v = v.replace(',', '')
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def parse_wide(path, sheet, codes, year):
    """Rows of an ONS wide table: code, sex, All ages, age groups. Older releases give the code on the male row
    only, so the last code seen is carried to the following female row. Returns long frame (code, sex, age, n)."""
    d = pd.read_excel(path, sheet_name=sheet, header=None, dtype=object)
    out, last = [], None
    allcode = re.compile(r'^[CD]\d\d(\.\d)?$')
    for _, row in d.iterrows():
        cells = row.tolist()
        strs = [(i, c.strip()) for i, c in enumerate(cells) if isinstance(c, str) and c.strip()]
        code_hit = next(((i, c) for i, c in strs if re.match(r'^[CD]\d\d', c)), None)   # codes and code ranges
        sex_hit = next(((i, c) for i, c in strs if c in ('M', 'F', 'Males', 'Females')), None)
        if code_hit is not None:
            last = code_hit[1]
        if sex_hit is None or last not in codes:
            continue
        if code_hit is None and sex_hit[1] in ('M', 'Males'):
            continue                                  # a male row always carries its code
        sex = 1 if sex_hit[1] in ('M', 'Males') else 2
        vals = [num(c) for c in cells[sex_hit[0] + 1:]]
        vals = [v for v in vals if not (isinstance(v, float) and np.isnan(v))]
        total = vals[0]
        for ages in (AGES20, AGES19):
            if len(vals) >= 1 + len(ages) and abs(sum(vals[1:1 + len(ages)]) - total) < 0.5:
                for a, v in zip(ages, vals[1:1 + len(ages)]):
                    out.append({'row': _, 'year': year, 'code': last, 'sex': sex, 'age': to_grid(a), 'n': v})
                break
        else:
            raise ValueError(f'age columns do not add up: {path.name} {sheet} {last} {sex}')
    o = pd.DataFrame(out).groupby(['row', 'year', 'code', 'sex', 'age'], as_index=False).n.sum()
    dup = o.groupby(['code', 'sex', 'age']).size()
    if (dup > 1).any():
        # some releases repeat an identical row (e.g. 2014 Table 1); accept exact repeats only
        g = o.groupby(['code', 'sex', 'age']).n.agg(['min', 'max'])
        bad = g[g['min'] != g['max']]
        if len(bad):
            raise ValueError(f'conflicting repeated rows: {path.name} {sheet} {sorted(set(bad.index.get_level_values(0)))}')
        o = o.drop_duplicates(['code', 'sex', 'age'])
    o = o.drop(columns='row')
    missing = {(c, s_) for c in codes for s_ in (1, 2)} - set(zip(o.code, o.sex))
    if missing:
        raise ValueError(f'missing rows: {path.name} {sheet} {sorted(missing)}')
    return o


def parse_2017(path):
    t = pd.read_excel(path, sheet_name='Table 1 Data', header=0, dtype=object)
    t = t[t['Registration Type'] == 'Number of Diagnoses']
    agecols = [c for c in t.columns if c not in ('Registration Type', 'ICD-10 code', 'Site description', 'Sex', 'All ages')]
    rows = []
    for _, r in t[t['ICD-10 code'].astype(str).str.strip().isin(CODES3)].iterrows():
        vals = [num(r[c]) for c in agecols]
        ages = AGES20 if len(agecols) == len(AGES20) else AGES19
        assert abs(sum(vals) - num(r['All ages'])) < 0.5
        for a, v in zip(ages, vals):
            rows.append({'year': 2017, 'code': str(r['ICD-10 code']).strip(), 'sex': 1 if r['Sex'] == 'Males' else 2,
                         'age': to_grid(a), 'n': v})
    return pd.DataFrame(rows).groupby(['year', 'code', 'sex', 'age'], as_index=False).n.sum()


def registrations():
    frames = []
    for y in YEARS:
        p = REG / (f'ons_{y}.xlsx' if y == 2017 else f'ons_{y}.xls')
        three = parse_2017(p) if y == 2017 else parse_wide(p, 'Table 1', set(CODES3), y)
        four = parse_wide(p, 'Table 8' if y in (2015, 2016) else 'Table 6', set(CODES4) | {'C71', 'D43'}, y)
        # consistency: fourth-character table totals for C71 and D43 equal the third-character table
        for c in ['C71', 'D43']:
            a = three[three.code == c].n.sum()
            b = four[four.code == c].n.sum()
            assert abs(a - b) < 0.5, (y, c, a, b)
        frames += [three, four[four.code.isin(CODES4)]]
    r = pd.concat(frames, ignore_index=True)
    assert set(r.year) == set(YEARS)
    return r


def deaths_ew():
    """WHO Mortality Database, England and Wales (4310): deaths by fourth-character code, sex and age."""
    col_age = {**{f'Deaths{i}': 0 for i in range(2, 7)}, **{f'Deaths{7 + i}': 5 * (i + 1) for i in range(16)},
               'Deaths23': 85, 'Deaths24': 90, 'Deaths25': 95}
    rows = []
    for part in range(1, 7):
        z = zipfile.ZipFile(RAW / 'who_mdb' / f'morticd10_part{part}.zip')
        d = pd.read_csv(z.open(z.namelist()[0]), dtype={'Cause': str, 'Country': str, 'List': str, 'Frmat': str}, low_memory=False)
        d = d[(d.Country == '4310') & d.Admin1.isna() & d.SubDiv.isna() & d.Year.between(2001, 2020) & d.Sex.isin([1, 2])]
        d = d[d.Cause.str[:3].isin(CODES3)]
        rows.append(d)
    d = pd.concat(rows, ignore_index=True)
    assert (d.List == '104').all(), d.List.unique()
    long = d.melt(id_vars=['Year', 'Sex', 'Cause'], value_vars=list(col_age), var_name='col', value_name='n')
    long['age'] = long.col.map(col_age).map(to_grid)
    long['n'] = long.n.fillna(0)
    out = []
    for code3 in CODES3:
        s = long[long.Cause.str[:3] == code3]
        out.append(s.groupby(['Year', 'Sex', 'age'], as_index=False).n.sum().assign(code=code3))
    for c4 in CODES4:
        s = long[long.Cause == c4.replace('.', '')]
        out.append(s.groupby(['Year', 'Sex', 'age'], as_index=False).n.sum().assign(code=c4))
    o = pd.concat(out, ignore_index=True).rename(columns={'Year': 'year', 'Sex': 'sex'})
    return o[['year', 'code', 'sex', 'age', 'n']]


def population():
    z = zipfile.ZipFile(REG / 'ons_pop_detailed_2001_2020.zip')
    name = 'MYEB1_detailed_population_estimates_series_UK_(2020_geog20).csv'
    p = pd.read_csv(z.open(name), low_memory=False)
    p = p[p.country.isin(['E', 'W'])]
    long = p.melt(id_vars=['country', 'sex', 'age'], value_vars=[f'population_{y}' for y in range(2001, 2021)],
                  var_name='year', value_name='pop')
    long['year'] = long.year.str[-4:].astype(int)
    long['age'] = long.age.map(to_grid)
    eng = long[long.country == 'E'].groupby(['year', 'sex', 'age'], as_index=False)['pop'].sum()
    ew = long.groupby(['year', 'sex', 'age'], as_index=False)['pop'].sum()
    return eng, ew


def group_counts(df):
    rows = []
    for g, codes in GROUPS.items():
        s = df[df.code.isin(codes)].groupby(['year', 'sex', 'age'], as_index=False).n.sum()
        rows.append(s.assign(group=g))
    return pd.concat(rows, ignore_index=True)


def asr(df, pop):
    w = esp_weights(GRID)
    w = w / w.sum()
    m = df.groupby(['measure', 'group', 'year', 'age'], as_index=False).n.sum().merge(
        pop.groupby(['popset', 'year', 'age'], as_index=False)['pop'].sum(), on=['year', 'age'], how='left')
    m = m[((m.measure == 'incidence') & (m.popset == 'E')) | ((m.measure == 'mortality') & (m.popset == 'EW'))]
    m['w'] = m.age.map(w)
    m['r'] = m.n / m['pop'] * m.w * 1e5
    out = m.groupby(['measure', 'group', 'year'], as_index=False).agg(asr=('r', 'sum'), n=('n', 'sum'))
    return out


def trend(df, pop, measure, group, years, step=None):
    s = df[(df.measure == measure) & (df.group == group) & df.year.isin(years)]
    ps = pop[pop.popset == ('E' if measure == 'incidence' else 'EW')]
    s = s.merge(ps, on=['year', 'sex', 'age'])
    s = s[s['pop'] > 0]
    X = pd.get_dummies(s[['age', 'sex']].astype(str), drop_first=True).astype(float)
    X['year_c'] = s.year - years[0]
    if step is not None:                       # level shift from the year of a coding change (sensitivity analysis)
        X['step'] = (s.year >= step).astype(float)
    X = sm.add_constant(X)
    fit = sm.GLM(s.n, X, family=sm.families.Poisson(), offset=np.log(s['pop'])).fit(scale='X2')
    b, se = fit.params['year_c'], fit.bse['year_c']
    out = {'measure': measure, 'group': group, 'years': f'{years[0]}-{years[-1]}', 'step_year': step, 'n_total': s.n.sum(),
           'aapc': 100 * (np.exp(b) - 1), 'lo': 100 * (np.exp(b - 1.96 * se) - 1), 'hi': 100 * (np.exp(b + 1.96 * se) - 1),
           'scale': fit.scale}
    if step is not None:
        bs_, ses_ = fit.params['step'], fit.bse['step']
        out.update({'step_pct': 100 * (np.exp(bs_) - 1), 'step_lo': 100 * (np.exp(bs_ - 1.96 * ses_) - 1),
                    'step_hi': 100 * (np.exp(bs_ + 1.96 * ses_) - 1)})
    return out


def trend_difference(df, pop, group, years):
    """Difference between the mortality and the incidence trend (ratio of annual rate ratios) from one Poisson model
    with measure-specific age-sex terms and a measure x year interaction; quasi-Poisson scale."""
    parts = []
    for m, ps in (('incidence', 'E'), ('mortality', 'EW')):
        s = df[(df.measure == m) & (df.group == group) & df.year.isin(years)].merge(pop[pop.popset == ps], on=['year', 'sex', 'age'])
        parts.append(s[s['pop'] > 0].assign(mort=float(m == 'mortality')))
    s = pd.concat(parts, ignore_index=True)
    X = pd.get_dummies((s.measure + '_' + s.age.astype(str) + '_' + s.sex.astype(str)).to_frame('cell'), drop_first=True).astype(float)
    X['year_c'] = s.year - years[0]
    X['year_mort'] = X.year_c * s.mort
    X = sm.add_constant(X)
    fit = sm.GLM(s.n, X, family=sm.families.Poisson(), offset=np.log(s['pop'])).fit(scale='X2')
    b, se = fit.params['year_mort'], fit.bse['year_mort']
    return {'group': group, 'years': f'{years[0]}-{years[-1]}', 'ratio_pct': 100 * (np.exp(b) - 1),
            'lo': 100 * (np.exp(b - 1.96 * se) - 1), 'hi': 100 * (np.exp(b + 1.96 * se) - 1)}


def main():
    reg = group_counts(registrations()).assign(measure='incidence')
    dth = group_counts(deaths_ew()).assign(measure='mortality')
    eng, ew = population()
    pop = pd.concat([eng.assign(popset='E'), ew.assign(popset='EW')], ignore_index=True)
    allc = pd.concat([reg, dth], ignore_index=True)
    allc.to_csv(TAB / 'registry_england_counts.csv', index=False)
    a = asr(allc, pop)
    a.to_csv(TAB / 'registry_england_asr.csv', index=False)
    tr = [trend(allc, pop, m, g, YEARS) for m in ['incidence', 'mortality'] for g in GROUPS]
    tr += [trend(allc, pop, 'mortality', g, list(range(2003, 2018))) for g in GROUPS]
    # sensitivity: level shift in 2014 (introduction of Iris in England and Wales; registry series also changed then)
    tr += [trend(allc, pop, m, g, YEARS, step=2014) for m in ['incidence', 'mortality'] for g in GROUPS]
    tr = pd.DataFrame(tr)
    tr.to_csv(TAB / 'registry_england_trends.csv', index=False)
    td = pd.DataFrame([trend_difference(allc, pop, g, YEARS) for g in GROUPS]
                      + [trend_difference(allc, pop, g, list(range(2010, 2018))) for g in GROUPS])
    td.to_csv(TAB / 'registry_england_trend_difference.csv', index=False)
    print(td.round(2).to_string(index=False))
    piv = a.pivot_table(index=['group', 'year'], columns='measure', values=['asr', 'n']).reset_index()
    piv.columns = ['_'.join([c for c in col if c]) for col in piv.columns]
    piv['mi_ratio'] = piv.asr_mortality / piv.asr_incidence
    sh = []
    for y in YEARS:
        g = a[a.year == y].set_index(['measure', 'group']).n
        sh.append({'year': y, 'c719_share_reg': g[('incidence', 'C71.9')] / g[('incidence', 'C71')],
                   'c719_share_dth': g[('mortality', 'C71.9')] / g[('mortality', 'C71')],
                   'd432_share_reg': g[('incidence', 'D43.2')] / g[('incidence', 'D43')],
                   'd432_share_dth': g[('mortality', 'D43.2')] / g[('mortality', 'D43')],
                   'bu_share_reg': g[('incidence', 'BU')] / g[('incidence', 'COMP')],
                   'bu_share_dth': g[('mortality', 'BU')] / g[('mortality', 'COMP')]})
    piv = piv.merge(pd.DataFrame(sh), on='year', how='left')
    piv.to_csv(TAB / 'registry_england_ratios.csv', index=False)
    pd.set_option('display.width', 200)
    print(a[a.year.isin([2009, 2017])].pivot_table(index='group', columns=['measure', 'year'], values='asr').round(2))
    print(tr.round(2).to_string(index=False))
    print(pd.DataFrame(sh).round(3).to_string(index=False))


if __name__ == '__main__':
    main()
