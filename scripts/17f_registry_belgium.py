"""Step 17f: cancer registry incidence by tumour behaviour, Belgium (deviations table, entry 22).

Question: how do deaths recorded with each group of CNS tumour codes compare with registrations of tumours with the
corresponding topography and behaviour in the same population? In particular, are deaths coded D43 (uncertain or
unknown behaviour) matched by registrations of CNS tumours of borderline behaviour, and did the documented changes of
coding system in 2012 and 2017 shift deaths without a corresponding change in registrations?

Inputs
- Belgian Cancer Registry. Primary brain and other central nervous system tumours in adults in Belgium, 2004-2020
  (data/raw/registry_belgium/bcr_cns_adults_2004_2020.pdf), appendix V: number of new diagnoses by primary location
  (ICD-O-3 topography C70, C71, C72, C75.1-C75.3), behaviour (malignant /3, borderline /1, benign /0) and incidence
  year, adults aged 20 years and older, both sexes. All CNS tumours have been notifiable since 2004; diagnoses
  without histological confirmation are included. The table is read from the PDF with pdftotext -layout.
- WHO Mortality Database, Belgium: deaths by underlying cause (fourth character), age 20 years and older, 2004-2020
  (analysis panel from step 9), and the population aged 20 years and older (Eurostat, as in the panel).

Correspondence between registry groups and ICD-10 mortality codes (pituitary, craniopharyngeal duct and pineal gland
are excluded, as in the mortality outcomes):
  malignant C70 + C71 + C72  -> C70-C72        borderline C70 -> D42        benign C70 -> D32
  borderline C71 + C72       -> D43            benign C71 + C72 -> D33
  borderline C71             -> D43.0-D43.2 (brain)

Outputs (output/tables)
- registry_belgium_counts.csv: registrations and deaths (age 20+) by group and year, with the population aged 20+
- registry_belgium_ratios.csv: ratio of deaths to registrations by group, pooled over the period and for its first and
  last three years (2004-2020 for malignant and borderline groups; 2010-2020 for groups containing benign tumours,
  whose registration was incomplete before 2010 according to the registry)
- registry_belgium_lymphoma.csv: malignant haematolymphoid tumours by location, 2004-2020 pooled (report table 2);
  the ratios file adds the malignant ratio with these tumours removed from registrations (pooled period only,
  because annual counts by location are not published)
- registry_belgium_survival.csv: five-year relative survival of borderline tumours by cohort (appendix VI)
- registry_belgium_trends.csv: average annual change over the same periods of crude rates at age 20+ (Poisson, log population
  offset, quasi-Poisson scale), with and without level shifts in 2012 and 2017, and the difference between the mortality and
  incidence trends (measure x year interaction); age-adjusted mortality trends for comparison with crude trends
"""
import re
import subprocess

import numpy as np
import pandas as pd
import statsmodels.api as sm

from config import RAW, OUTPUT
from lib.data import load_panel

SRC = RAW / 'registry_belgium' / 'bcr_cns_adults_2004_2020.pdf'
TAB = OUTPUT / 'tables'
YEARS = list(range(2004, 2021))
STEPS = [2012, 2017]                          # documented coding changes: Iris in Flanders and Brussels (2012) and in
                                              # Wallonia, replacing STYX (2017) (Eurostat metadata; step terms in step 12)
TITLE = ('Belgium: Number of new diagnoses and age-standardised incidence of primary brain and other CNS tumours in '
         'adults by primary location, behaviour and incidence year (2004-2020)')
LOCS = {'Tumours of the meninges': 'C70', 'Tumours of the brain': 'C71',
        'Tumours of the spinal cord, cranial nerves and other': 'C72',
        'Tumours of the pituitary and pineal glands and': 'C75', 'All primary brain and other CNS tumours': 'ALL'}
BEH = {'Malignant behaviour': 'mal', 'Borderline behaviour': 'bord', 'Benign behaviour': 'ben'}
GROUPS = {'MAL': [('C70', 'mal'), ('C71', 'mal'), ('C72', 'mal')], 'D42': [('C70', 'bord')], 'D32': [('C70', 'ben')],
          'D43': [('C71', 'bord'), ('C72', 'bord')], 'D43_brain': [('C71', 'bord')], 'D33': [('C71', 'ben'), ('C72', 'ben')]}
GROUPS['BU'] = GROUPS['D42'] + GROUPS['D32'] + GROUPS['D43'] + GROUPS['D33']
GROUPS['COMP'] = GROUPS['MAL'] + GROUPS['BU']
BENIGN = {'D32', 'D33', 'BU', 'COMP'}          # contain benign tumours, registered completely from 2010 (BCR report)
DEATH_COLS = {'MAL': ['MAL'], 'D42': ['D42'], 'D32': ['D32'], 'D43': ['D43'], 'D43_brain': ['D43_2', 'D43_01'],
              'D33': ['D33'], 'BU': ['D32', 'D33', 'D42', 'D43'], 'COMP': ['MAL', 'D32', 'D33', 'D42', 'D43']}


def read_table():
    """Parse appendix V (males and females) from the layout text; columns located from the year header."""
    txt = subprocess.run(['pdftotext', '-layout', str(SRC), '-'], capture_output=True, text=True, check=True).stdout
    lines = txt.splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip() == TITLE and 'Males and females' in lines[i + 1])
    end = next(i for i in range(start + 1, len(lines)) if lines[i].strip().startswith('Belgium:'))
    block = lines[start:end]
    hdr = next(l for l in block if re.search(r'2004\s+2005', l))
    assert [int(x) for x in re.findall(r'20\d\d', hdr)[:len(YEARS)]] == YEARS

    def values(l):
        """First 17 cells after the label: counts separated by two or more spaces (thousands separated by one space)."""
        toks = re.split(r'\s{2,}', l.strip())
        while toks and not re.fullmatch(r'\d{1,3}(?: \d{3})*|-', toks[0]):
            toks = toks[1:]
        toks = toks[:len(YEARS)]
        assert len(toks) == len(YEARS) and all(re.fullmatch(r'\d{1,3}(?: \d{3})*|-', x) for x in toks), l
        return [0 if x == '-' else int(x.replace(' ', '')) for x in toks]

    rows, loc = [], None
    for i, l in enumerate(block):
        s = l.strip()
        lab = re.split(r'\s{2,}', s)[0] if s else ''
        indent = len(l) - len(l.lstrip())
        if lab in LOCS and indent <= 2:
            loc = LOCS[lab]
            num_line = l if re.search(r'\s{2,}\d', l.strip()) else block[i + 1]
            rows.append({'loc': loc, 'beh': 'all', **dict(zip(YEARS, values(num_line)))})
        elif lab in BEH and loc and indent <= 5:
            rows.append({'loc': loc, 'beh': BEH[lab], **dict(zip(YEARS, values(l)))})
    t = pd.DataFrame(rows).set_index(['loc', 'beh'])
    assert len(t) == 20, len(t)
    for lc in ['C70', 'C71', 'C72', 'C75', 'ALL']:                       # behaviours add up to each location total
        assert (t.loc[(lc, 'all')] == t.loc[[(lc, b) for b in BEH.values()]].sum()).all(), lc
    for b in list(BEH.values()) + ['all']:                               # locations add up to the overall total
        assert (t.loc[('ALL', b)] == t.loc[[(lc, b) for lc in ['C70', 'C71', 'C72', 'C75']]].sum()).all(), b
    return t


def read_survival():
    """Five-year relative survival by cohort for tumours of the brain and all CNS tumours of borderline behaviour
    (appendix VI, males and females)."""
    txt = subprocess.run(['pdftotext', '-layout', str(SRC), '-'], capture_output=True, text=True, check=True).stdout
    lines = txt.splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip().startswith('Belgium: 5-year relative survival trends')
                 and 'Males and females' in lines[i + 1])
    out, loc = [], None
    for l in lines[start:start + 80]:
        lab = re.split(r'\s{2,}', l.strip())[0]
        if lab in ('Tumours of the brain', 'All primary brain and other CNS tumours'):
            loc = lab
        elif lab == 'Borderline behaviour' and loc and len(l) - len(l.lstrip()) <= 5:
            toks = re.split(r'\s{2,}', l.strip())[1:]
            out.append({'location': loc, 'n_2004_2009': int(toks[0].replace(' ', '')), 'n_2010_2015': int(toks[1].replace(' ', '')),
                        'n_2016_2020': int(toks[2].replace(' ', '')), 'rs5_2004_2009': float(toks[3].replace(',', '.')),
                        'rs5_2010_2015': float(toks[4].replace(',', '.')), 'rs5_2016_2020': float(toks[5].replace(',', '.'))})
            loc = None
        if len(out) == 2:
            break
    assert len(out) == 2, out
    return pd.DataFrame(out)


def read_lymphoma():
    """Malignant haematolymphoid tumours (mostly lymphomas) and all malignant tumours by location, 2004-2020 pooled
    (report table 2, first column). Deaths from these tumours are coded to lymphoma (C81-C96), not to C70-C72."""
    txt = subprocess.run(['pdftotext', '-layout', str(SRC), '-'], capture_output=True, text=True, check=True).stdout
    lines = txt.splitlines()
    start = next(i for i, l in enumerate(lines) if re.match(r'\s*Table 2\s+All primary brain and other CNS tumours in adults', l))
    heads = {'Meninges (C70)': 'C70', 'Brain (C71)': 'C71', 'Spinal cord, cranial nerves and other parts of CNS (C72)': 'C72'}
    out, loc = {}, None
    for l in lines[start:start + 60]:
        st = l.strip()
        if st in heads:
            loc = heads[st]
        elif '(C75' in st:
            loc = None
        elif loc and re.match(r'(All tumours|Haematolymphoid tumours)\s{2,}', st):
            lab, n = re.split(r'\s{2,}', st)[:2]
            assert (loc, lab) not in out
            out[(loc, lab)] = int(n.replace(',', ''))
    assert len(out) == 6, out
    return pd.DataFrame([{'loc': k, 'all_malignant': out[(k, 'All tumours')], 'haematolymphoid': out[(k, 'Haematolymphoid tumours')]}
                         for k in heads.values()])


def build():
    t = read_table()
    p = load_panel()
    p = p[(p.geo == 'BE') & (p.age >= 20) & p.year.isin(YEARS)]
    dth = p.groupby('year')[sorted({c for v in DEATH_COLS.values() for c in v})].sum()
    pop = p.groupby('year')['pop'].sum()
    rows = []
    for g in GROUPS:
        for y in YEARS:
            rows.append({'group': g, 'year': y, 'registrations': int(sum(t.loc[k, y] for k in GROUPS[g])),
                         'deaths': int(dth.loc[y, DEATH_COLS[g]].sum()), 'pop20': float(pop.loc[y])})
    return pd.DataFrame(rows), p


def years_for(g):
    return list(range(2010, 2021)) if g in BENIGN else YEARS


def ratios(c):
    out = []
    for g, s in c.groupby('group', sort=False):
        yy0 = years_for(g)
        for lab, yy in [('all', yy0), ('first3', yy0[:3]), ('last3', yy0[-3:])]:
            q = s[s.year.isin(yy)]
            d, r = q.deaths.sum(), q.registrations.sum()
            # 95% CI of the ratio of two Poisson counts (log scale)
            se = np.sqrt(1 / d + 1 / r) if d > 0 and r > 0 else np.nan
            out.append({'group': g, 'period': lab, 'years': f'{yy[0]}-{yy[-1]}', 'deaths': d, 'registrations': r, 'ratio': d / r,
                        'lo': d / r * np.exp(-1.96 * se), 'hi': d / r * np.exp(1.96 * se)})
    return pd.DataFrame(out)


def fit_trend(s, col, step=False):
    X = pd.DataFrame({'year_c': s.year - s.year.min()}, index=s.index)
    if step:
        for y in STEPS:
            X[f'step{y}'] = (s.year >= y).astype(float)
    X = sm.add_constant(X)
    f = sm.GLM(s[col], X, family=sm.families.Poisson(), offset=np.log(s.pop20)).fit(scale='X2')
    r = {'aapc': 100 * (np.exp(f.params.year_c) - 1), 'lo': 100 * (np.exp(f.params.year_c - 1.96 * f.bse.year_c) - 1),
         'hi': 100 * (np.exp(f.params.year_c + 1.96 * f.bse.year_c) - 1)}
    if step:
        for y in STEPS:
            b_, se = f.params[f'step{y}'], f.bse[f'step{y}']
            r.update({f'step{y}_pct': 100 * (np.exp(b_) - 1), f'step{y}_lo': 100 * (np.exp(b_ - 1.96 * se) - 1),
                      f'step{y}_hi': 100 * (np.exp(b_ + 1.96 * se) - 1)})
    return r


def trend_difference(s):
    a = s.assign(n=s.registrations, mort=0.0)
    b = s.assign(n=s.deaths, mort=1.0)
    d = pd.concat([a, b], ignore_index=True)
    X = pd.DataFrame({'mort': d.mort, 'year_c': d.year - d.year.min()})
    X['year_mort'] = X.year_c * X.mort
    X = sm.add_constant(X)
    f = sm.GLM(d.n, X, family=sm.families.Poisson(), offset=np.log(d.pop20)).fit(scale='X2')
    b_, se = f.params.year_mort, f.bse.year_mort
    return {'diff_pct': 100 * (np.exp(b_) - 1), 'diff_lo': 100 * (np.exp(b_ - 1.96 * se) - 1),
            'diff_hi': 100 * (np.exp(b_ + 1.96 * se) - 1)}


def asr_mortality_trend(p, cols, years):
    """Average annual change in mortality adjusted for age (Poisson with five-year age terms, ages 20+)."""
    q = p[p.year.isin(years)].groupby(['year', 'age'], as_index=False)[cols + ['pop']].sum()
    q['n'] = q[cols].sum(axis=1)
    X = pd.get_dummies(q.age.astype(str), drop_first=True).astype(float)
    X['year_c'] = q.year - years[0]
    X = sm.add_constant(X)
    f = sm.GLM(q.n, X, family=sm.families.Poisson(), offset=np.log(q['pop'])).fit(scale='X2')
    return {'aapc_age_adj': 100 * (np.exp(f.params.year_c) - 1),
            'lo_age_adj': 100 * (np.exp(f.params.year_c - 1.96 * f.bse.year_c) - 1),
            'hi_age_adj': 100 * (np.exp(f.params.year_c + 1.96 * f.bse.year_c) - 1)}


def main():
    c, p = build()
    c.to_csv(TAB / 'registry_belgium_counts.csv', index=False)
    r = ratios(c)
    ly = read_lymphoma()
    t = read_table()
    for _, x in ly.iterrows():                 # table 2 totals agree with appendix V, 2004-2020
        assert x.all_malignant == t.loc[(x['loc'], 'mal'), YEARS].sum(), x
    q = r[(r.group == 'MAL') & (r.period == 'all')].iloc[0]
    d, n = q.deaths, q.registrations - ly.haematolymphoid.sum()
    se = np.sqrt(1 / d + 1 / n)
    r = pd.concat([r, pd.DataFrame([{'group': 'MAL_excl_lymphoma', 'period': 'all', 'years': q.years, 'deaths': d,
                                     'registrations': n, 'ratio': d / n, 'lo': d / n * np.exp(-1.96 * se),
                                     'hi': d / n * np.exp(1.96 * se)}])], ignore_index=True)
    r.to_csv(TAB / 'registry_belgium_ratios.csv', index=False)
    ly.to_csv(TAB / 'registry_belgium_lymphoma.csv', index=False)
    sv = read_survival()
    sv.to_csv(TAB / 'registry_belgium_survival.csv', index=False)
    print(sv.to_string(index=False))
    tr = []
    for g, s in c.groupby('group', sort=False):
        s = s[s.year.isin(years_for(g))].reset_index(drop=True)
        row = {'group': g, 'years': f'{years_for(g)[0]}-{years_for(g)[-1]}'}
        for m, col in [('inc', 'registrations'), ('mort', 'deaths')]:
            row.update({f'{m}_{k}': v for k, v in fit_trend(s, col).items()})
            row.update({f'{m}_s_{k}': v for k, v in fit_trend(s, col, step=True).items()})
        row.update(trend_difference(s))
        row.update(asr_mortality_trend(p, DEATH_COLS[g], years_for(g)))
        tr.append(row)
    tr = pd.DataFrame(tr)
    tr.to_csv(TAB / 'registry_belgium_trends.csv', index=False)
    pd.set_option('display.width', 250)
    print(c.pivot_table(index='year', columns='group', values=['registrations', 'deaths']).to_string())
    print(r.round(2).to_string(index=False))
    print(tr.round(1).to_string(index=False))


if __name__ == '__main__':
    main()
