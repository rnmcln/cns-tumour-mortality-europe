# Source data and how to reconstruct them

The analysis uses only published aggregate data. None of the source files are included in this repository, because they are distributed by their publishers under their own terms. This file lists each source, the files the scripts expect, how to obtain them and the dates on which the files used for the published analysis were downloaded.

All files go under `data/raw/`, in the subfolders given below. The download scripts in `scripts/` create these folders and write the files with the expected names; `./run_all.sh --fetch` runs all of them. `data/SHA256SUMS_raw_inputs.txt` contains the SHA-256 checksums of the files used for the published analysis (78 files). From the repository root:

```
sha256sum -c data/SHA256SUMS_raw_inputs.txt
```

The scripts do not require the checksums to match, except `02c_fetch_registries.py`, which reports any mismatch for the registry files. Where a publisher has updated a file since the access date, results may differ from the published ones.

## Overview

| Source | Folder | Script | Downloaded (UTC) | Can change after download |
|---|---|---|---|---|
| WHO Mortality Database, ICD-10 raw data files (release of February 2026) | `who_mdb/` | `01_fetch_who_mdb.py` | 29 October 2026 | yes, with each WHO release |
| Eurostat: deaths, populations, EUROPOP2023, life expectancy | `eurostat/` | `02_fetch_eurostat.py` | 29 October 2026 | yes (live datasets) |
| Eurostat: EUROPOP2025 | `eurostat/` | `02b_fetch_europop2025.py` | 30 October 2026 | yes (live dataset) |
| UN World Population Prospects 2024 | `unwpp/` | `02b_fetch_unwpp_metadata.py` | 30 October 2026 | no (fixed revision) |
| Eurostat reference metadata on causes of death, with annexes | `metadata/` | `02b_fetch_unwpp_metadata.py` | 30 October 2026 | yes |
| Office for National Statistics (ONS), cancer registrations, England, 2009–2017, and population estimates | `registry_england/` | `02c_fetch_registries.py` | 30 October 2026 | unlikely (archived releases) |
| Belgian Cancer Registry, report on CNS tumours in adults, 2004–2020 | `registry_belgium/` | `02c_fetch_registries.py` | 30 October 2026 | unlikely (published report) |
| Coding changes in the UK and Germany (authors' compilation) | `metadata/additional_coding_changes.csv` | included in this repository | | |

Each download script writes an `ACCESSED.txt` file with the date and time of download.

## 1. WHO Mortality Database

Deaths by country, year, sex, five-year age group and underlying cause (ICD-10), and the populations reported to WHO.

- Download page: https://www.who.int/data/data-collection-tools/who-mortality-database
- Files (folder `who_mdb/`): `morticd10_part1.zip` to `morticd10_part6.zip`, `morticd10_add.zip` (restricted-format files, used for Germany from 2021 and Norway from 2017), `mort_pop.zip`, `mort_country_codes.zip`, `mort_availability.zip`, `mort_notes.zip` and the documentation file `mort_documentation71f9e29d-7e3f-41e6-aafc-c4c1775c7aa3.zip`.
- `01_fetch_who_mdb.py` downloads the files from `https://cdn.who.int/media/docs/default-source/world-health-data-platform/mortality-raw-data/`. The file URLs carry a version suffix (listed in the script) that changes with each WHO release. If a URL no longer resolves, download the current files from the download page, save them under the names above and run `python3 01_fetch_who_mdb.py --no-download`. A later release will include revised and additional data, and results will then differ from the published ones.
- The script reads cause list 103, 104 and 10M records for the study countries and writes `data/interim/who_deaths.csv` and `data/interim/who_population.csv`.

## 2. Eurostat

Downloaded with the `eurostat` Python package (https://pypi.org/project/eurostat/), which queries the Eurostat dissemination API. Datasets (folder `eurostat/`, one gzip-compressed csv file each):

| Dataset | Content | Filter used |
|---|---|---|
| `hlth_cd_aro` | causes of death by country of residence | ICD-10 C70-C72, C, TOTAL, R, D00-D48; residents; numbers |
| `demo_pjan` | population on 1 January by age and sex | numbers |
| `proj_23np` | EUROPOP2023 population projections, all variants | persons; Serbia not covered |
| `proj_25np` | EUROPOP2025 population projections (sensitivity analysis) | persons; Serbia and the UK not covered |
| `demo_mlexpec` | life expectancy by age and sex | years |

Data browser: https://ec.europa.eu/eurostat/databrowser/ (for example https://ec.europa.eu/eurostat/databrowser/view/hlth_cd_aro). Eurostat revises these datasets, so a new download is not expected to reproduce the checksums. The analysis uses the WHO data for deaths; the Eurostat death counts are used for checks and sensitivity analyses.

## 3. UN World Population Prospects 2024

Population by five-year age group and sex, medium variant, used for projections for Serbia and the UK.

- File: `WPP2024_PopulationByAge5GroupSex_Medium.csv.gz` from https://population.un.org/wpp/downloads (CSV format, "Population by 5-year age groups").
- `02b_fetch_unwpp_metadata.py` downloads this file, keeps the rows for Serbia and the UK as `unwpp/wpp2024_srb_gbr_age5.csv.gz` and deletes the full file.

## 4. Eurostat reference metadata on causes of death

National reference metadata (SIMS) for the cause-of-death statistics, used to identify coding systems (for example Iris or MMDS) and the years in which they were introduced.

- Pages: `https://ec.europa.eu/eurostat/cache/metadata/EN/hlth_cdeath_simscd_<cc>.htm`, where `<cc>` is the lower-case Eurostat country code; saved as `metadata/hlth_cdeath_simscd_<cc>.htm`.
- Annexes linked from these pages: `https://ec.europa.eu/eurostat/cache/metadata/Annexes/<file>.xlsx`; saved in `metadata/annexes/`.
- The overall metadata page `metadata/hlth_cdeath_sims.htm` (https://ec.europa.eu/eurostat/cache/metadata/EN/hlth_cdeath_sims.htm) is listed in the checksums for reference; no script reads it.
- Downloaded by `02b_fetch_unwpp_metadata.py`. Eurostat updates these pages, so the content may differ from that used for the published analysis.

The coding changes for the UK and Germany that are not documented in these pages were compiled from national sources and are included as `metadata/additional_coding_changes.csv`; each row gives its source.

## 5. Office for National Statistics, England

Cancer registrations in England by site and behaviour for diagnosis years 2009 to 2017, and population estimates.

- Dataset page: https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/cancerregistrationstatisticscancerregistrationstatisticsengland
- One file per diagnosis year, saved as `registry_england/ons_<year>.xls` (2017: `ons_2017.xlsx`). The source file names on ons.gov.uk are listed in `02c_fetch_registries.py`.
- Population estimates, mid-2001 to mid-2020 detailed time series (MYEB1), saved as `registry_england/ons_pop_detailed_2001_2020.zip`, from https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland
- Registrations for 2018 onwards are published by NHS England (National Disease Registration Service) and were not used.

## 6. Belgian Cancer Registry

Belgian Cancer Registry. Primary brain and other central nervous system tumours in adults in Belgium 2004–2020. Brussels: Belgian Cancer Registry, 2023 (D/2023/11.846/2).

- Landing page: https://kankerregister.org/en/publicaties/primary-brain-and-other-central-nervous-system-tumours-adults-belgium-2004-2020
- File: https://kankerregister.org/sites/default/files/2024/bcrpub_neuroonco_en_0.pdf, saved as `registry_belgium/bcr_cns_adults_2004_2020.pdf`.
- `17f_registry_belgium.py` reads appendix V (new diagnoses by location, behaviour and incidence year) and appendix VI (five-year relative survival) from the text layer with `pdftotext -layout`.

## Sources assessed and not used

Netherlands: the public NKR Cijfers tool of the Netherlands Comprehensive Cancer Organisation (https://nkr-cijfers.iknl.nl) reports incidence of CNS tumours by topography but not by behaviour, so it could not be used for a comparison by behaviour. Counts by behaviour would require a data request to IKNL.

## Derived data

The pipeline writes intermediate files to `data/interim/` and `data/processed/`, and results to `output/`. These folders are not part of the repository and are recreated by `run_all.sh`.
