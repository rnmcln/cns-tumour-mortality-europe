"""Step 8: detailed WHO extract for the main analysis (protocol v1.0).

Reads the raw WHO ICD-10 files (parts 1-6 and the restricted-format file for Germany 2021+ and
Norway 2017+) and keeps, for the 32 countries (national level only), the rows needed for the
outcomes and indicators of the protocol:
  AAA (all causes); every C code (C00-C97); every D code (D00-D48; D32, D33, D42, D43 for the outcomes, the others for the additional analyses); every R code.
All original columns are kept (age columns Deaths1-Deaths26 unchanged, suppressed cells left blank),
so that age harmonisation and gap-filling are done transparently in step 9.
Output: data/interim/who_detail.csv.gz (one row per country, year, list, cause, sex).
"""
import zipfile

import pandas as pd

from config import RAW, INTERIM, WHO_CODES

SRC = RAW / 'who_mdb'
STR = {'Admin1': str, 'SubDiv': str, 'List': str, 'Cause': str, 'Frmat': str, 'IM_Frmat': str}
COLS = ['Country', 'Year', 'List', 'Cause', 'Sex', 'Frmat'] + [f'Deaths{i}' for i in range(1, 27)]


def keep_cause(c):
    c = str(c)
    if c == 'AAA':
        return True
    if c[0] == 'C' and c[1:3].isdigit() and int(c[1:3]) <= 97:
        return True
    if c[0] == 'D' and c[1:3].isdigit() and int(c[1:3]) <= 48:   # D00-D48; non-CNS codes for the negative control
        return True
    return c[0] == 'R'


def national(df):
    keep = {str(k) for k in WHO_CODES}
    return df[df.Country.astype(str).isin(keep) & df.Admin1.isna() & df.SubDiv.isna()]


def main():
    frames = []
    for part in range(1, 7):
        with zipfile.ZipFile(SRC / f'morticd10_part{part}.zip') as z:
            with z.open(z.namelist()[0]) as fh:
                for chunk in pd.read_csv(fh, dtype=str, chunksize=300_000):
                    chunk = national(chunk)
                    mask = chunk.Cause.map(keep_cause).astype(bool).to_numpy()
                    frames.append(chunk.loc[mask, COLS].assign(src_file=f'part{part}'))
    with zipfile.ZipFile(SRC / 'morticd10_add.zip') as z:
        for n in z.namelist():
            sheet = 'mortICD_DEU' if 'Germany' in n else 'mortICD_NOR'
            x = pd.read_excel(z.open(n), sheet_name=sheet, dtype=STR)
            x = national(x.astype({'Country': str}))
            mask = x.Cause.map(keep_cause).astype(bool).to_numpy()
            frames.append(x.loc[mask, COLS].assign(src_file=f'add_{sheet[-3:]}'))
    d = pd.concat(frames, ignore_index=True)
    d['Country'] = d.Country.astype(int)
    d['geo'] = d.Country.map(WHO_CODES)
    d['Year'] = d.Year.astype(int)
    num = [f'Deaths{i}' for i in range(1, 27)]
    d[num] = d[num].apply(pd.to_numeric, errors='coerce')
    # a country-year must come from one file only (the restricted file supersedes nothing in parts 1-6)
    nfile = d.groupby(['geo', 'Year']).src_file.nunique()
    assert nfile.max() == 1, nfile[nfile > 1]
    d = d.sort_values(['geo', 'Year', 'Cause', 'Sex']).reset_index(drop=True)
    d.to_csv(INTERIM / 'who_detail.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
    print('rows:', len(d), '| country-years:', d.groupby(['geo', 'Year']).ngroups)


if __name__ == '__main__':
    main()
