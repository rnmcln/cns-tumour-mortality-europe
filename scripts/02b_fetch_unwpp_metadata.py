"""Step 2b (needs internet): UN World Population Prospects 2024 (Serbia, United Kingdom) and the
Eurostat country metadata on cause-of-death coding (pages and xlsx annexes).

Outputs: data/raw/unwpp/wpp2024_srb_gbr_age5.csv.gz; data/raw/metadata/*.htm and annexes/*.xlsx.
"""
import datetime as dt
import re
import urllib.request

import pandas as pd

from config import RAW, COUNTRIES

STAMP = dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%MZ') + '\n'

# UN WPP 2024, population by five-year age group and sex, medium variant (1950-2100, 1 July, thousands)
WPP = ('https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/'
       'WPP2024_PopulationByAge5GroupSex_Medium.csv.gz')
out = RAW / 'unwpp'
out.mkdir(parents=True, exist_ok=True)
full = out / 'WPP2024_PopulationByAge5GroupSex_Medium.csv.gz'
urllib.request.urlretrieve(WPP, full)
d = pd.read_csv(full, low_memory=False)
d[d.ISO3_code.isin(['SRB', 'GBR'])].to_csv(out / 'wpp2024_srb_gbr_age5.csv.gz', index=False,
                                             compression={'method': 'gzip', 'mtime': 0})
full.unlink()                      # only the two-country extract is kept
(out / 'ACCESSED.txt').write_text(STAMP)

# Eurostat country metadata for causes of death (SIMS), with xlsx annexes
meta = RAW / 'metadata'
(meta / 'annexes').mkdir(parents=True, exist_ok=True)
BASE = 'https://ec.europa.eu/eurostat/cache/metadata/'
for cc in COUNTRIES:
    name = f'hlth_cdeath_simscd_{cc.lower()}.htm'
    urllib.request.urlretrieve(f'{BASE}EN/{name}', meta / name)
    page = (meta / name).read_text(encoding='utf-8', errors='ignore')
    for a in sorted(set(re.findall(r'\.\./Annexes/(hlth_cdeath_simscd_[^"]+\.xlsx)', page))):
        urllib.request.urlretrieve(f'{BASE}Annexes/{a}', meta / 'annexes' / a)
(meta / 'ACCESSED.txt').write_text(STAMP)
print('done')
