# Tumour behaviour coding, population ageing and central nervous system tumour mortality in 26 European countries, 2003–22: a population-based study

Analysis code and documentation for the study of the same title.

**Authors**

Aaron Lawson McLean (ORCID 0000-0001-5528-6905), Julian Kahr, Anne Neumeister, Christian Senft

1. Department of Neurosurgery, Jena University Hospital, Friedrich Schiller University Jena, Jena, Germany
2. Comprehensive Cancer Centre Central Germany, Jena Campus, Jena, Germany

All authors have both affiliations.

**Correspondence:** Aaron Lawson McLean, Neuro-Oncology and Health Intelligence Lab, Department of Neurosurgery, Jena University Hospital, Friedrich Schiller University Jena, Am Klinikum 1, 07747 Jena, Germany. aaron.lawsonmclean@med.uni-jena.de

**Protocol:** version 1.0, deposited on Zenodo before the main analysis, [doi:10.5281/zenodo.23797629](https://doi.org/10.5281/zenodo.23797629). A copy is in `protocol/`.

## Summary

European statistics on central nervous system (CNS) tumour mortality count malignant codes (ICD-10 C70–C72) only. Deaths certified as CNS tumours of benign or uncertain behaviour (D32–D33, D42–D43) are excluded, and their share depends on certification practice. The study assessed how this definition affects trends and country comparisons, and the contribution of population ageing.

Using WHO vital registration data for 26 countries, age-standardised mortality (2013 European Standard Population) and country ranks were compared for malignant codes and for all CNS tumour codes. The change in malignant deaths between 2003–05 and 2020–22 was decomposed into population size, age structure, change in all-CNS mortality and a non-malignant component. Projection models were evaluated out of sample, and registry data from England and Belgium were compared with the mortality data.

Age-standardised mortality from all CNS tumours fell from 8.92 to 8.32 per 100,000, whereas malignant mortality was unchanged (6.69 and 6.64). Country ranks on the two definitions correlated moderately (Spearman 0.59). Malignant deaths increased from 28,767 to 34,547 per year; population ageing and growth contributed +6,278. A decline in deaths with non-malignant codes, mainly brain tumours of uncertain behaviour, added +2,622, of which about 1,690 (95% interval 860 to 2,240) was matched by higher malignant mortality. Projected malignant deaths in 2040 ranged from 39,365 to 42,114 across scenarios.

All numbers above are produced by the code in this repository; the tables that contain them are listed under [Outputs](#outputs).

## Repository structure

```
protocol/                      protocol v1.0 (md, pdf); country sets; published list of deviations
scripts/                       analysis pipeline (Python and R), numbered in run order
scripts/config.py              paths, country lists, ICD-10 code groups, age groups, seeds
scripts/lib/                   shared functions (data access, decomposition, rate methods, projection models)
data/README.md                 how to reconstruct the source data
data/SHA256SUMS_raw_inputs.txt checksums of the source files used for the published analysis
data/raw/metadata/additional_coding_changes.csv
                               coding changes in the UK and Germany compiled from national sources (authors' compilation)
run_all.sh                     runs the full pipeline
requirements.txt               Python packages
R_packages.txt                 R packages and system tools
```

No source data are included, because they are published by third parties under their own terms. `data/README.md` lists every source file, where to obtain it and how the scripts use it.

## Reproducing the analysis

1. Install Python 3.11 or later with the packages in `requirements.txt`, and R 4.3 with the packages in `R_packages.txt` (INLA and BAPC are not on CRAN; see the links in that file). Step 17f needs `pdftotext` from poppler-utils.
2. Obtain the source data, either automatically:

   ```
   ./run_all.sh --fetch
   ```

   which downloads all files into `data/raw/` and then runs the analysis, or manually as described in `data/README.md`, followed by

   ```
   ./run_all.sh
   ```

3. Compare the downloaded files with `data/SHA256SUMS_raw_inputs.txt` (`sha256sum -c data/SHA256SUMS_raw_inputs.txt` from the repository root). Files that are updated by their publishers (see `data/README.md`) may differ from the versions used for the published analysis.

A full run from the source files took about 21 minutes in a test on a single core, most of it in the R models (steps 12 and 16).

### Pipeline steps

| Step | Script | Content |
|---|---|---|
| 1 | `01_fetch_who_mdb.py` | WHO Mortality Database: download and long-format extract for the study countries |
| 2 | `02_fetch_eurostat.py`, `02b_fetch_europop2025.py`, `02b_fetch_unwpp_metadata.py`, `02c_fetch_registries.py` | Eurostat deaths, populations, projections and life expectancy; EUROPOP2025; UN World Population Prospects 2024; Eurostat coding metadata; ONS and Belgian Cancer Registry files |
| 3 | `03_harmonise.py` | harmonised deaths and populations by country, year, sex and five-year age group |
| 4–5 | `04_pilot.py`, `05_pilot_figures.py` | feasibility pilot run before protocol deposition (not part of `run_all.sh`) |
| 6 | `06_coding_metadata.py` | coding systems and change years from the Eurostat metadata |
| 7 | `07_country_sets.py` | country sets for the primary and secondary analyses |
| 8–9 | `08_extract_who_detail.py`, `09_build_panel.py` | ICD-10 subcode detail; analysis panel |
| 10 | `10_coding_changes.py` | documented coding changes per country |
| 11 | `11_decomposition.py` | decomposition of the change in deaths (Das Gupta / Shapley), with bootstrap intervals |
| 12 | `12_trend_model.R` | hierarchical negative binomial trend models (glmmTMB), with validation fits ending in 2010, 2013 and 2016 |
| 13 | `13_rates_joinpoint.py` | age-standardised rates and joinpoint trends |
| 14 | `14_population_projections.py` | population projections to 2040 |
| 15–16 | `15_projections.py`, `16_bapc.R` | projection models, rolling-origin evaluation and projections to 2040; Bayesian age-period-cohort models (BAPC/INLA) |
| 17 | `17_secondary.py`, `17b`–`17f` | secondary and sensitivity analyses; coding evidence; forecast checks; comparisons between definitions; registry comparisons for England and Belgium |
| 18 | `18_figures.py` | figures |

### Outputs

The pipeline writes tables to `output/tables/`, model objects to `output/models/`, figures to `output/figures/` and the software versions used to `output/software_versions.txt`. The main tables are:

| Table | Content |
|---|---|
| `asr_by_period_setL.csv` | age-standardised rates by period, primary country set |
| `comparison_country_ranks.csv`, `comparison_rank_summary.csv`, `comparison_trend_direction.csv` | country ranks and trends on malignant and all CNS tumour codes |
| `decomp_primary.csv`, `decomp_primary_bootstrap.csv` | primary decomposition |
| `coding_evidence_*.csv` | non-malignant component and the estimate of transfer to malignant codes |
| `trend_*.csv` | trend models and out-of-sample evaluation |
| `projections_aggregate.csv`, `projections_country.csv`, `backtest_*.csv` | projections to 2040 and projection backtests |
| `registry_england_*.csv`, `registry_belgium_*.csv` | registry comparisons |

### Reproducibility notes

- Random seeds are fixed in `scripts/config.py` and in the R scripts. Python and BLAS are run single-threaded (`run_all.sh`).
- The R model fits are not bitwise reproducible between runs. Estimates from the glmmTMB trend models (step 12) differ from run to run in about the sixth significant digit, and the BAPC/INLA predictions (step 16) by a median of about one part in a million and at most about 0.5% in strata with few deaths. In a full re-run from the source files, all numbers in the main text of the paper were reproduced exactly; some simulated interval bounds and one summary of the Bayesian age-period-cohort backtest in the appendix differed in the last digit.
- The analysis departs from the protocol in the respects listed in `protocol/deviations_published.csv`.

## Licence

Code: MIT licence (`LICENSE`). The source data remain subject to the terms of their publishers (see `data/README.md`).

## Citation

Please cite the study and this repository; `CITATION.cff` contains the repository metadata. The protocol can be cited as: Lawson McLean A. Demographic change, death certification and recorded mortality from central nervous system tumours in Europe: study protocol for a decomposition of change since 2003 and projections to 2040. Version 1.0. Zenodo, 2026. doi:10.5281/zenodo.23797629.
