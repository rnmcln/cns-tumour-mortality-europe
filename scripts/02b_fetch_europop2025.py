"""Step 2b: EUROPOP2025 population projections (Eurostat proj_25np), for a sensitivity analysis of the
projections (protocol/deviations_published.csv, entry 19). Writes data/raw/eurostat/proj_25np.csv.gz and
ACCESSED_proj_25np.txt. EUROPOP2023 remains the population basis of the main projections (protocol)."""
import datetime as dt

import eurostat

from config import RAW, COUNTRIES

OUT = RAW / 'eurostat'
geo = [g for g in COUNTRIES if g not in ('RS', 'UK')]
df = eurostat.get_data_df('proj_25np', filter_pars={'geo': geo, 'unit': 'PER'})
df.to_csv(OUT / 'proj_25np.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
(OUT / 'ACCESSED_proj_25np.txt').write_text(dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%MZ') + '\n')
print('proj_25np', df.shape, sorted(df.projection.unique()), df['geo\\TIME_PERIOD'].nunique())
