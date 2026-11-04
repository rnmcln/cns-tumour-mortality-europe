"""Step 1: WHO Mortality Database (ICD-10 files, population, documentation).

Downloads the raw files (needs internet), then writes a filtered long-format extract for the
31 study countries: data/interim/who_deaths.csv and data/interim/who_population.csv.
Run with --no-download to rebuild the extract from files already in data/raw/who_mdb.
"""
import datetime as dt
import sys
import urllib.request
import zipfile

import pandas as pd

from config import RAW, INTERIM, WHO_CODES, AGE_LOWER

SRC = RAW / 'who_mdb'
BASE = 'https://cdn.who.int/media/docs/default-source/world-health-data-platform/mortality-raw-data/'
FILES = {  # file name -> version suffix as published on the WHO download page
    'mort_availability.zip': '23261e11_33', 'mort_country_codes.zip': '800faac2_5',
    'mort_documentation71f9e29d-7e3f-41e6-aafc-c4c1775c7aa3.zip': '40cce9be_42',
    'mort_notes.zip': '2020c859_19', 'mort_pop.zip': '937039fc_26', 'morticd10_add.zip': '47111c9f_4',
    'morticd10_part1.zip': 'e2a4f93a_19', 'morticd10_part2.zip': '6e55000b_6',
    'morticd10_part3.zip': '9f1111a2_14', 'morticd10_part4.zip': '259c5c23_30',
    'morticd10_part5.zip': 'ad970d0b_34', 'morticd10_part6.zip': 'ec801a61_4'}

# Deaths1 = all ages; Deaths2 = 0; Deaths3-6 = 1..4; Deaths7 = 5-9; ...; Deaths23 = 85-89 (or 85+
# in formats 01/02); Deaths24 = 90-94; Deaths25 = 95+; Deaths26 = unknown age.
GROUP_COLS = {0: ['Deaths2', 'Deaths3', 'Deaths4', 'Deaths5', 'Deaths6']}
for i, a in enumerate(AGE_LOWER[1:-1]):
    GROUP_COLS[a] = [f'Deaths{7 + i}']
GROUP_COLS[85] = ['Deaths23', 'Deaths24', 'Deaths25']
STR = {'Admin1': str, 'SubDiv': str, 'List': str, 'Cause': str, 'Frmat': str, 'IM_Frmat': str}


def cause_group(c):
    """Map an ICD-10 code from lists 103/104/10M to the analysis groups (None = not needed)."""
    if c == 'AAA':
        return 'ALL'
    c3 = c[:3]
    if c3 in ('C70', 'C71', 'C72'):
        return 'C70-C72'
    if c3[0] == 'C' and c3 <= 'C97':
        return 'C_OTHER'          # malignant neoplasms other than C70-C72
    if c3[0] == 'R':
        return 'R00-R99'
    if c3 in ('D32', 'D33', 'D42', 'D43'):
        return 'D32-D33_D42-D43'
    return None


def download():
    SRC.mkdir(parents=True, exist_ok=True)
    for name, ver in FILES.items():
        urllib.request.urlretrieve(f'{BASE}{name}?sfvrsn={ver}', SRC / name)
    (SRC / 'ACCESSED.txt').write_text(dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%MZ') + '\n')


def read_deaths():
    frames = []
    keep = {str(c) for c in WHO_CODES}
    for part in range(1, 7):                  # filter each part on reading to limit memory use
        with zipfile.ZipFile(SRC / f'morticd10_part{part}.zip') as z:
            with z.open(z.namelist()[0]) as fh:
                for chunk in pd.read_csv(fh, dtype=str, chunksize=200_000):
                    frames.append(chunk[chunk.Country.isin(keep)])
    with zipfile.ZipFile(SRC / 'morticd10_add.zip') as z:
        for n in z.namelist():   # Germany 2021+, Norway 2017+ (restricted formats)
            frames.append(pd.read_excel(z.open(n), sheet_name=f'mortICD_{"DEU" if "Germany" in n else "NOR"}',
                                        dtype=STR))
    d = pd.concat(frames, ignore_index=True)
    d['Country'] = d['Country'].astype(int)
    d['Year'] = d['Year'].astype(int)
    d = d[d.Country.isin(WHO_CODES) & d.Admin1.isna() & d.SubDiv.isna()].copy()
    d['group'] = d.Cause.map(cause_group)
    d = d[d.group.notna()]
    d['geo'] = d.Country.map(WHO_CODES)
    num = [f'Deaths{i}' for i in range(1, 27)]
    d[num] = d[num].apply(pd.to_numeric, errors='coerce')
    rows = []
    for (geo, year, sex, frmat, grp), s in d.groupby(['geo', 'Year', 'Sex', 'Frmat', 'group']):
        tot = s[num].sum(min_count=1)
        rec = {'geo': geo, 'year': year, 'sex': int(sex), 'frmat': frmat, 'cause': grp,
               'deaths_all_ages': tot['Deaths1'], 'deaths_unknown_age': tot['Deaths26'],
               'n_codes': len(s)}
        for a, cols in GROUP_COLS.items():
            rec[f'a{a}'] = tot[cols].sum(min_count=1)
        rows.append(rec)
    return pd.DataFrame(rows)


def read_population():
    with zipfile.ZipFile(SRC / 'mort_pop.zip') as z:
        p = pd.read_csv(z.open(z.namelist()[0]), dtype=STR)
    p['Country'] = p['Country'].astype(int)
    p = p[p.Country.isin(WHO_CODES) & p.Admin1.isna() & p.SubDiv.isna()].copy()
    p['geo'] = p.Country.map(WHO_CODES)
    num = [f'Pop{i}' for i in range(1, 27)]
    p[num] = p[num].apply(pd.to_numeric, errors='coerce')
    out = p[['geo', 'Year', 'Sex', 'Frmat']].rename(columns={'Year': 'year', 'Sex': 'sex', 'Frmat': 'frmat'})
    out['pop_all_ages'] = p['Pop1'].values
    for a, cols in GROUP_COLS.items():
        out[f'a{a}'] = p[[c.replace('Deaths', 'Pop') for c in cols]].sum(axis=1, min_count=1).values
    return out


if __name__ == '__main__':
    if '--no-download' not in sys.argv:
        download()
    INTERIM.mkdir(parents=True, exist_ok=True)
    deaths = read_deaths()
    deaths.to_csv(INTERIM / 'who_deaths.csv', index=False)
    read_population().to_csv(INTERIM / 'who_population.csv', index=False)
    print('WHO deaths records:', len(deaths))
