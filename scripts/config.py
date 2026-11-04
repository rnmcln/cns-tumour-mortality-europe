"""Shared configuration for the trends and projections study."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data' / 'raw'
INTERIM = ROOT / 'data' / 'interim'
PROCESSED = ROOT / 'data' / 'processed'
PILOT = ROOT / 'pilot'
OUTPUT = ROOT / 'output'

# Same 31 countries as the 2023 burden study (EU-27, Iceland, Norway, Serbia, Switzerland)
COUNTRIES = ['AT', 'BE', 'BG', 'CZ', 'HR', 'CY', 'DK', 'EE', 'FI', 'FR', 'DE', 'EL', 'HU', 'IS', 'IE',
             'IT', 'LV', 'LT', 'LU', 'MT', 'NL', 'NO', 'PL', 'PT', 'RO', 'RS', 'SI', 'SK', 'ES', 'SE', 'CH']

# WHO Mortality Database country codes
WHO_CODES = {4010: 'AT', 4020: 'BE', 4030: 'BG', 4038: 'HR', 3080: 'CY', 4045: 'CZ', 4050: 'DK',
             4055: 'EE', 4070: 'FI', 4080: 'FR', 4085: 'DE', 4140: 'EL', 4150: 'HU', 4160: 'IS',
             4170: 'IE', 4180: 'IT', 4186: 'LV', 4188: 'LT', 4190: 'LU', 4200: 'MT', 4210: 'NL',
             4220: 'NO', 4230: 'PL', 4240: 'PT', 4270: 'RO', 4273: 'RS', 4274: 'SK', 4276: 'SI',
             4280: 'ES', 4290: 'SE', 4300: 'CH', 4308: 'UK'}

# Study countries plus the United Kingdom (extension; WHO data and ONS/UN WPP populations)
COUNTRIES_EXT = COUNTRIES + ['UK']

# Harmonised age groups: lower bounds of 0-4, 5-9, ..., 80-84, 85+ (common to WHO formats 00-02)
AGE_LOWER = list(range(0, 90, 5))
AGE_LABELS = [f'{a}-{a + 4}' for a in AGE_LOWER[:-1]] + ['85+']

# 2013 European Standard Population, collapsed to the harmonised groups (85+ = 85-89 + 90-94 + 95+)
ESP2013 = [5000, 5500, 5500, 5500, 6000, 6000, 6500, 7000, 7000, 7000, 7000, 6500, 6000, 5500,
           5000, 4000, 2500, 1500 + 800 + 200]

CAUSE_BRAIN = 'C70-C72'
FIRST_YEAR_WHO, FIRST_YEAR_EUROSTAT, LAST_YEAR = 2000, 2011, 2023
BASE_PERIOD = (2019, 2023)     # rates held constant in the demographic scenario
PROJECTION_YEARS = (2024, 2040)
SEED = 20261029

# ---- Main analysis (protocol v1.0) -------------------------------------------------------------
# UN M49 sub-regions of Europe; Cyprus (M49: Western Asia) is assigned to Southern Europe
REGION = {**{g: 'Northern Europe' for g in ['DK', 'EE', 'FI', 'IS', 'IE', 'LV', 'LT', 'NO', 'SE', 'UK']},
          **{g: 'Western Europe' for g in ['AT', 'BE', 'FR', 'DE', 'LU', 'NL', 'CH']},
          **{g: 'Southern Europe' for g in ['HR', 'CY', 'EL', 'IT', 'MT', 'PT', 'RS', 'SI', 'ES']},
          **{g: 'Eastern Europe' for g in ['BG', 'CZ', 'HU', 'PL', 'RO', 'SK']}}
NAMES = {'AT': 'Austria', 'BE': 'Belgium', 'BG': 'Bulgaria', 'HR': 'Croatia', 'CY': 'Cyprus', 'CZ': 'Czechia',
         'DK': 'Denmark', 'EE': 'Estonia', 'FI': 'Finland', 'FR': 'France', 'DE': 'Germany', 'EL': 'Greece',
         'HU': 'Hungary', 'IS': 'Iceland', 'IE': 'Ireland', 'IT': 'Italy', 'LV': 'Latvia', 'LT': 'Lithuania',
         'LU': 'Luxembourg', 'MT': 'Malta', 'NL': 'Netherlands', 'NO': 'Norway', 'PL': 'Poland', 'PT': 'Portugal',
         'RO': 'Romania', 'RS': 'Serbia', 'SK': 'Slovakia', 'SI': 'Slovenia', 'ES': 'Spain', 'SE': 'Sweden',
         'CH': 'Switzerland', 'UK': 'United Kingdom'}
# 2013 European Standard Population, five-year groups to 95+ (lower bound -> weight)
ESP2013_FULL = dict(zip(list(range(0, 100, 5)), [5000, 5500, 5500, 5500, 6000, 6000, 6500, 7000, 7000, 7000,
                                                 7000, 6500, 6000, 5500, 5000, 4000, 2500, 1500, 800, 200]))
# elementary death groups in the analysis panel
GROUPS = ['MAL', 'MAL_UNSPEC', 'MAL_3CHAR', 'D32', 'D33', 'D42', 'D43', 'CALL', 'R', 'ALL',
          'D43_2', 'D43_01', 'D43_OTH', 'C71_9', 'C79_3', 'C80', 'R90_0', 'SELLAR', 'DALL']
WINDOW = (2000, 2023)          # years whose age formats fix each country's age grid
BASE_L = (2003, 2005)
END_L = (2020, 2022)
BASE_W = (2011, 2013)
END_W = (2021, 2023)
N_BOOT = 2000
