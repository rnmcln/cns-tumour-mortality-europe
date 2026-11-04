"""Step 3: harmonise WHO and Eurostat deaths and populations to one long format.

Output (data/processed):
  deaths_long.csv      source, geo, year, sex (M/F; Eurostat also T = both sexes as published), cause,
                       age (lower bound), deaths
  population_long.csv  source, geo, year, sex (M/F; Eurostat also T), age, population (Eurostat: 1 January; WHO: mid-year)
Age groups (lower bounds): WHO 0-4, ..., 80-84, 85+ (finest grouping common to WHO formats 00-02);
Eurostat 0-4, ..., 90-94, 95+. Analyses that combine sources collapse Eurostat to 85+.
"""
import pandas as pd

from config import RAW, INTERIM, PROCESSED, AGE_LOWER, COUNTRIES

PROCESSED.mkdir(parents=True, exist_ok=True)
SEX = {1: 'M', 2: 'F', 'M': 'M', 'F': 'F'}


def eurostat_age(code):
    """Map Eurostat five-year age codes to harmonised lower bounds (None = aggregate code)."""
    if code in ('Y_LT1', 'Y1-4'):
        return 0
    if code == 'Y_GE95':
        return 95
    if code.startswith('Y') and '-' in code:
        lo, hi = code[1:].split('-')
        if int(hi) - int(lo) == 4:
            return int(lo)
    return None


def eurostat_deaths():
    d = pd.read_csv(RAW / 'eurostat' / 'hlth_cd_aro.csv.gz')
    d = d.rename(columns={'geo\\TIME_PERIOD': 'geo'})
    d = d[d.sex.isin(['M', 'F', 'T']) & d.geo.isin(COUNTRIES)]
    years = [c for c in d.columns if c.isdigit()]
    long = d.melt(id_vars=['geo', 'sex', 'age', 'icd10'], value_vars=years, var_name='year',
                  value_name='deaths').dropna(subset=['deaths'])
    long['year'] = long.year.astype(int)
    long['age_h'] = long.age.map(eurostat_age)
    long = long[long.age_h.notna()]
    out = long.groupby(['geo', 'year', 'sex', 'icd10', 'age_h'], as_index=False).deaths.sum()
    out = out.rename(columns={'icd10': 'cause', 'age_h': 'age'})
    out['cause'] = out.cause.replace({'TOTAL': 'ALL', 'C': 'C00-C97', 'R': 'R00-R99'})
    out['source'] = 'eurostat'
    return out


def who_deaths():
    d = pd.read_csv(INTERIM / 'who_deaths.csv', dtype={'frmat': str})
    d = d[d.sex.isin([1, 2])]
    agec = [f'a{a}' for a in AGE_LOWER]
    long = d.melt(id_vars=['geo', 'year', 'sex', 'cause', 'frmat'], value_vars=agec, var_name='age',
                  value_name='deaths')
    long['age'] = long.age.str[1:].astype(int)
    long['deaths'] = long.deaths.fillna(0)
    long['sex'] = long.sex.map(SEX)
    # all malignant neoplasms = C70-C72 + other C codes
    c = long[long.cause.isin(['C70-C72', 'C_OTHER'])].groupby(['geo', 'year', 'sex', 'age', 'frmat'],
                                                               as_index=False).deaths.sum()
    c['cause'] = 'C00-C97'
    long = pd.concat([long[long.cause != 'C_OTHER'], c], ignore_index=True)
    long['source'] = 'who'
    return long


def eurostat_population():
    p = pd.read_csv(RAW / 'eurostat' / 'demo_pjan.csv.gz').rename(columns={'geo\\TIME_PERIOD': 'geo'})
    p = p[p.sex.isin(['M', 'F', 'T']) & p.geo.isin(COUNTRIES)]
    years = [c for c in p.columns if c.isdigit() and int(c) >= 1999]

    def single(a):
        if a == 'Y_LT1':
            return 0
        if a == 'Y_OPEN':
            return 100
        if a.startswith('Y') and a[1:].isdigit():
            return int(a[1:])
        return None
    p['a1'] = p.age.map(single)
    p = p[p.a1.notna()]
    long = p.melt(id_vars=['geo', 'sex', 'a1'], value_vars=years, var_name='year',
                  value_name='population').dropna(subset=['population'])
    long['year'] = long.year.astype(int)
    long['age'] = (long.a1 // 5 * 5).clip(upper=95).astype(int)
    out = long.groupby(['geo', 'year', 'sex', 'age'], as_index=False).population.sum()
    out['source'] = 'eurostat'
    return out


def who_population():
    p = pd.read_csv(INTERIM / 'who_population.csv', dtype={'frmat': str})
    p = p[p.sex.isin([1, 2])]
    agec = [f'a{a}' for a in AGE_LOWER]
    long = p.melt(id_vars=['geo', 'year', 'sex', 'frmat'], value_vars=agec, var_name='age',
                  value_name='population')
    long['age'] = long.age.str[1:].astype(int)
    long['sex'] = long.sex.map(SEX)
    long['source'] = 'who'
    return long


if __name__ == '__main__':
    deaths = pd.concat([eurostat_deaths(), who_deaths()], ignore_index=True)
    deaths.to_csv(PROCESSED / 'deaths_long.csv', index=False)
    pop = pd.concat([eurostat_population(), who_population()], ignore_index=True)
    pop.to_csv(PROCESSED / 'population_long.csv', index=False)
    print(deaths.groupby(['source', 'cause']).deaths.sum().round(0).to_string())
