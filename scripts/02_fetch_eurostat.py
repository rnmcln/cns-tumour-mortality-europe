"""Step 2: Eurostat extracts (needs internet). Writes data/raw/eurostat/*.csv.gz and ACCESSED.txt."""
import datetime as dt
import sys

import eurostat

from config import RAW, COUNTRIES

OUT = RAW / 'eurostat'
OUT.mkdir(parents=True, exist_ok=True)
GEO = COUNTRIES

SKIP = set(sys.argv[1:])          # e.g. "python 02_fetch_eurostat.py hlth_cd_aro demo_pjan" skips those


def get(code, **flt):
    if code in SKIP:
        return
    df = eurostat.get_data_df(code, filter_pars=flt)
    df.to_csv(OUT / f'{code}.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
    print(code, df.shape)


# Deaths 2011-2023: brain/CNS, all malignant neoplasms, all causes; by sex and five-year age group
get('hlth_cd_aro', geo=GEO, icd10=['C70-C72', 'C', 'TOTAL', 'R', 'D00-D48'], resid='TOT_IN', unit='NR')
# Population on 1 January 2000-2024 by single year of age and sex
get('demo_pjan', geo=GEO, unit='NR')
# EUROPOP2023 population projections (all variants) by single year of age and sex
# EUROPOP2023 does not cover Serbia; Serbia needs UN World Population Prospects (see protocol)
get('proj_23np', geo=[g for g in GEO if g != 'RS'], unit='PER')
# Life expectancy by age and sex (for the reference life tables)
get('demo_mlexpec', geo=GEO, unit='YR')
(OUT / 'ACCESSED.txt').write_text(dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%MZ') + '\n')
