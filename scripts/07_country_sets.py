"""Step 7: declare the analytic country sets (protocol v1.0, section 4) from data availability only.

Rules
- Set L (long series, both outcomes): WHO deaths reported in every year of 2003-2005 and 2020-2022,
  and a denominator for each of these years (Eurostat 1 January populations; UK: WHO-reported ONS
  mid-year populations, 2021-2022 interpolated).
- Set W (wide, malignant outcome): deaths in every year of 2011-2013 and 2021-2023 (Eurostat; UK from WHO).
- Complete-case subset of Set L: countries with no gap-filled or restricted-file country-years in the two
  periods (excludes Germany 2021-2022 and Norway 2020-2022, which come from the WHO restricted file).
- Country-level reporting threshold: mean of at least 100 C70-C72 deaths per year in 2017-2019.
Output: protocol/country_sets.csv
"""
import pandas as pd

from config import INTERIM, PROCESSED, ROOT, COUNTRIES_EXT

w = pd.read_csv(INTERIM / 'who_deaths.csv', dtype={'frmat': str})
P = pd.read_csv(PROCESSED / 'population_long.csv')
D = pd.read_csv(PROCESSED / 'deaths_long.csv', dtype={'frmat': str})

who_c = w[w.cause == 'C70-C72'].groupby('geo').year.apply(set).to_dict()
eu_pop = P[(P.source == 'eurostat') & (P.sex == 'T')].groupby('geo').year.apply(set).to_dict()
uk_pop = set(P[(P.source == 'who') & (P.geo == 'UK')].year) | {2021, 2022}   # 2021-2022 interpolated
eu_c = D[(D.source == 'eurostat') & (D.sex == 'T') & (D.cause == 'C70-C72')].groupby('geo').year.apply(set).to_dict()
restricted = {('DE', y) for y in (2021, 2022)} | {('NO', y) for y in range(2017, 2025)}

L_YEARS = list(range(2003, 2006)) + list(range(2020, 2023))
W_YEARS = list(range(2011, 2014)) + list(range(2021, 2024))
rows = []
for g in COUNTRIES_EXT:
    pop = uk_pop if g == 'UK' else eu_pop.get(g, set())
    in_L = all(y in who_c.get(g, set()) and y in pop for y in L_YEARS)
    src_w = who_c.get(g, set()) if g == 'UK' else eu_c.get(g, set())
    in_W = all(y in src_w and y in pop for y in W_YEARS)
    cc = in_L and not any((g, y) in restricted for y in L_YEARS)
    if g == 'UK':
        base = w[(w.geo == 'UK') & (w.cause == 'C70-C72') & w.year.between(2017, 2019) & w.sex.isin([1, 2])]
        mean_c = base.deaths_all_ages.sum() / 3
    else:
        base = D[(D.source == 'eurostat') & (D.sex == 'T') & (D.cause == 'C70-C72') & (D.geo == g) & D.year.between(2017, 2019)]
        mean_c = base.deaths.sum() / 3
    rows.append({'geo': g, 'set_L': in_L, 'set_L_complete_case': cc, 'set_W': in_W,
                 'mean_C70_C72_deaths_2017_2019': round(mean_c, 1), 'country_level_reporting': mean_c >= 100})
out = pd.DataFrame(rows)
out.to_csv(ROOT / 'protocol' / 'country_sets.csv', index=False)
print(out.to_string())
print('Set L:', int(out.set_L.sum()), '| complete case:', int(out.set_L_complete_case.sum()), '| Set W:', int(out.set_W.sum()),
      '| country-level reporting:', int(out.country_level_reporting.sum()))
