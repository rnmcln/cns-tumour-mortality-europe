---
title: "Demographic change, death certification and recorded mortality from central nervous system tumours in Europe: study protocol for a decomposition of change since 2003 and projections to 2040"
version: "1.0"
date: "30 October 2026"
author: "Aaron Lawson McLean"
---

# Demographic change, death certification and recorded mortality from central nervous system tumours in Europe: study protocol

**Version:** 1.0, 30 October 2026

**Author:** Aaron Lawson McLean (ORCID 0000-0001-5528-6905), Department of Neurosurgery, Jena University Hospital, Friedrich Schiller University Jena, Am Klinikum 1, 07747 Jena, Germany. aaron.lawsonmclean@med.uni-jena.de

**Study type:** secondary analysis of published, aggregate vital-registration and population data. No individual-level data are used; ethical approval is not required.

## 1. Background

Recorded deaths from malignant neoplasms of the brain and other parts of the central nervous system (CNS) have increased in Europe while age-standardised mortality rates have been stable or have fallen. Analyses based on Global Burden of Disease (GBD) estimates attribute increases in high-income regions mainly to population ageing and growth, and project them forward (Kuang et al. 2025; Yuan et al. 2026). Long-term projections of cancer deaths by site are also published by the International Agency for Research on Cancer and the European Cancer Information System.

These analyses rest on modelled estimates. GBD redistributes deaths assigned to ill-defined codes and borrows strength across locations. It therefore cannot show how much of a recorded trend reflects changes in death certification and coding rather than changes in mortality.

In European vital-registration data, deaths from CNS tumours are assigned to malignant codes (C70–C72) or to codes for benign and uncertain-behaviour tumours (D32–D33, D42–D43), and the balance between the two changes over time. Some changes are abrupt and coincide with the introduction of automated coding systems; others are gradual. Malignant-code mortality trends therefore combine changes in disease, diagnosis and certification. Deaths from benign and uncertain-behaviour CNS tumours, which include meningiomas and tumours certified without confirmation of their behaviour, are not reported separately in European cancer-mortality statistics. To our knowledge, their contribution to recorded CNS tumour mortality has not been described for European countries individually.

## 2. Objectives

**Primary objective.** To partition the change in recorded CNS tumour deaths between 2003–2005 and 2020–2022 into the contributions of population size, population age-sex structure, coding change and residual change in age-specific mortality, by country, by region and for Europe.

**Secondary objectives.**

1. To describe the timing of these contributions (chained annual decomposition) and to estimate level shifts at documented and data-identified coding changes.
2. To estimate trends in age-specific mortality of the composite outcome, 2003–2019, and the difference between observed and expected deaths in 2020–2023.
3. To describe deaths from benign and uncertain-behaviour CNS tumours by age, sex, country and period, and between-country differences in their ratio to malignant CNS deaths.
4. To project deaths to 2040 under constant age-specific rates, and to 2035 under trend models, with the projection model selected by out-of-sample evaluation.
5. To describe CNS tumours as a share of all deaths from malignant neoplasms (C00–C97) over time.

## 3. Outcomes and measures

| Outcome | Definition (underlying cause of death, ICD-10) |
|---|---|
| Malignant outcome | C70–C72 |
| Composite outcome | C70–C72, D32–D33, D42–D43 |
| Composite outcome, narrow variant | C70–C72, D33, D43 |
| Benign and uncertain-behaviour CNS tumours | D32–D33, D42–D43 |

- **Measures:** numbers of deaths (primary); age-standardised mortality rates per 100,000 (2013 European Standard Population) for description; the share of deaths at ages 75+ and 85+.
- **Coding component:** the rate component of the malignant outcome minus the rate component of the composite outcome (section 5.1).
- **Certification indicators:** share of unspecified subcodes (C70.9, C71.9, C72.9) among malignant deaths; share of all deaths at ages 75+ assigned to R00–R99; ratio of benign and uncertain-behaviour to malignant CNS tumour deaths.

Secondary CNS neoplasms (C79.3) are excluded; deaths of patients with brain metastases are coded to the primary site. Codes are matched by their first three characters, so that three- and four-character reporting are treated alike.

## 4. Data

| Source | Content | Use |
|---|---|---|
| WHO Mortality Database, ICD-10 files, February 2026 release, including the file with restricted country formats | deaths by country, year, sex, age group and underlying cause | both outcomes; backbone of the long series |
| Eurostat `hlth_cd_aro` (extracted 29 October 2026) | deaths, European shortlist, 2011–2023 | malignant outcome for the wide country set; validation |
| Eurostat `demo_pjan` (extracted 29 October 2026) | population on 1 January by single year of age and sex | denominators (person-years as the mean of 1 January of years t and t+1) |
| Office for National Statistics mid-year estimates, as reported to the WHO Mortality Database | UK population | UK denominators |
| Eurostat `proj_23np` (EUROPOP2023) | projected population, baseline and five variants | projections for the 30 study countries other than Serbia |
| United Nations World Population Prospects 2024, medium variant | population by five-year age group and sex | projections for Serbia and the UK; sensitivity analysis for UK denominators |
| Eurostat country metadata on causes of death (`hlth_cdeath_simscd_<country>`, pages and annexes, retrieved 30 October 2026) | ICD-10 updates and automated coding systems by year, 2011 onwards | documented coding changes (Appendix A) |

**Countries.** The 27 EU member states, Iceland, Norway, Serbia and Switzerland, and the United Kingdom. Results are shown with and without the UK.

**Age groups.** Five-year groups to 95+ where the source provides them, otherwise to 85+. Where a country-year is reported to 85+ only, both periods of a comparison for that country use 85+. Norway's restricted-format data from 2017 combine ages 0–14; ages 0–14 are combined for Norway in all periods.

**Country sets** (fixed from data availability before any analysis of the primary outcome; Appendix B):

- **Set L (long series, both outcomes):** WHO deaths for every year of 2003–2005 and 2020–2022, with denominators. 26 countries: Austria, Belgium, Czechia, Denmark, Estonia, Finland, France, Germany, Hungary, Iceland, Italy, Latvia, Lithuania, Luxembourg, Malta, the Netherlands, Norway, Poland, Romania, Serbia, Slovakia, Slovenia, Spain, Sweden, Switzerland and the UK. The 25 countries other than the UK accounted for 89.7% of C70–C72 deaths in the 31 non-UK countries in 2019.
- **Complete-case subset of Set L:** Set L without Germany and Norway, whose 2020–2022 data come partly from restricted-format files (24 countries).
- **Set W (wide, malignant outcome):** all 32 countries, 2011–2013 against 2021–2023 (Eurostat; UK from the WHO database).
- **Country-level reporting:** countries with at least 100 C70–C72 deaths per year in 2017–2019 (28 countries). Iceland, Luxembourg, Malta and Cyprus contribute to regional and European results only.

**Gap-filling.**

- Germany 2021–2022 (cells with fewer than 3 deaths suppressed): the difference between the published all-ages total and the sum of age-specific cells is allocated to the suppressed age-sex cells in proportion to their share in 2018–2020.
- UK 2021–2022 denominators (missing in the WHO population file): geometric interpolation by age and sex between 2020 and 2023.
- No other gap-filling is performed.

**Data rules.** Only complete country-years are used. One ICD list per country-year. Eurostat both-sexes series are never summed with the sex-specific series.

## 5. Analysis

### 5.1 Decomposition (primary analysis)

For each country and outcome, the change in mean annual deaths between 2003–2005 and 2020–2022 is decomposed into population size, age-sex structure and age-sex-specific rates. The method is the symmetric decomposition of Das Gupta (1993), equivalent to the Shapley value over all orderings of the three factors. Country components are summed to UN M49 regions (Northern, Western, Southern and Eastern Europe) and to the European total, so no country-composition term arises.

The **coding component** is the difference between the rate components of the malignant and composite outcomes. The **residual rate component** is the rate component of the composite outcome.

A chained annual decomposition, accumulating year-to-year contributions, shows when each component accrued. Set W provides the same decomposition for the malignant outcome, 2011–2013 to 2021–2023.

Uncertainty is expressed in two ways:

- **Stochastic:** Poisson bootstrap of death counts, 2,000 replicates, fixed seed.
- **Systematic:** the range of results across the pre-specified sensitivity analyses (section 5.6), reported as the main expression of uncertainty.

### 5.2 Coding changes

Documented changes are the years in which the automated-coding status or the family of the coding system changed (Appendix A, 2011 onwards; software version updates alone are not counted).

Before 2011, a year t is a data-identified candidate change in a country with at least 100 C70–C72 deaths per year if both of the following hold:

- the absolute difference between the mean log ratio of benign and uncertain to malignant CNS deaths in years t to t+2 and in years t−3 to t−1 exceeds 0.3;
- the malignant and the benign or uncertain death rates change in opposite directions between these two windows.

At each documented or candidate change, a level shift is estimated jointly for the malignant and the benign or uncertain series (step terms in the model of section 5.3). Where an opposite-direction shift is estimated, the composite outcome is the primary outcome for that country in all analyses, and the malignant outcome is reported as secondary.

### 5.3 Trends

Age-sex-specific deaths, 2003–2019, are modelled with a hierarchical negative binomial model:

- person-year offset;
- fixed effects for age group and sex;
- country random intercepts and random slopes for calendar year, by sex and broad age group (0–39, 40–64, 65–74, 75–84, 85+);
- step terms at the coding changes from section 5.2.

Reported are the pooled slopes, partially pooled country slopes and the between-country variance; no country-specific significance tests are performed. Expected deaths for 2020–2023 are predicted from the fitted model and compared with observed deaths. Joinpoint regression is used only descriptively, with at most two joinpoints and a minimum segment of four years.

### 5.4 Projections

Candidate models:

1. constant age-sex-specific rates of 2015–2019 (benchmark);
2. age-specific log-linear trends with drift damped by 0, 25, 50 and 75% in successive five-year periods (Møller et al. 2003);
3. an age-drift model without cohort terms;
4. Nordpred, for countries with at least four observed five-year periods;
5. Bayesian age-period-cohort models (Riebler and Held 2017).

**Model selection.** Rolling-origin evaluation is used, with origins 2009, 2011, 2013 and 2015 and horizons of 5 and 10 years, for both outcomes and for countries with data from 2000 or earlier. Metrics are the absolute percentage error of total deaths and of deaths at 70+, the signed error and prediction-interval coverage; target years 2020–2022 are reported separately. A trend model replaces the constant-rate benchmark only if its error is lower at both horizons.

**Populations.**

- EUROPOP2023 baseline, with the lower-mortality and the three migration variants as range.
- Serbia and the UK: UN WPP 2024 age-sex-specific growth factors applied to the last observed population.
- Projected populations for years already observed are replaced by observed values.

**Horizons and uncertainty.** Horizons are 2040 for the constant-rate scenario and 2035 for trend models. Poisson, parameter and scenario uncertainty are reported separately.

### 5.5 Secondary analyses

- Deaths from benign and uncertain-behaviour CNS tumours by age group, sex, country and period.
- The ratio of these deaths to malignant CNS deaths by country and period.
- The share of CNS tumours among C00–C97 deaths over time, with the decomposition of C00–C97 deaths as reference.
- Certification indicators (section 3) by age group and period.

### 5.6 Pre-specified sensitivity analyses

1. Narrow composite outcome (D33 and D43 only).
2. Age groups truncated at 85+ for all countries.
3. End period 2017–2019 (before the COVID-19 pandemic).
4. Trends excluding 2020–2022.
5. Constant-rate projection with 2021–2023 rates.
6. WHO and Eurostat counts compared where both exist.
7. WHO population denominators instead of Eurostat, where available.
8. Exclusion of countries with documented or data-identified coding changes.
9. Complete-case subset of Set L.
10. Results with and without the UK.
11. Population projection variants; UN WPP denominators for the UK.

### 5.7 Software

Python 3 (pandas, NumPy, statsmodels, matplotlib) and R 4.3 (glmmTMB 1.1.8 for the hierarchical models; BAPC 0.0.37 with R-INLA 26.8.7 for the Bayesian age-period-cohort models). Installation of all packages was tested before this protocol was deposited. Fixed random seeds; package versions recorded with the code.

## 6. Interpretation

The residual rate component combines incidence, case fatality, diagnosis, certification and competing mortality, and is not interpreted as a change in risk. Case fatality of glioblastoma, the most frequent malignant CNS tumour in adults, is high, so mortality largely follows incidence and diagnosis. Improvements in survival are therefore expected to change mortality little.

## 7. Feasibility work before this protocol

Before this protocol was finalised, a feasibility pilot was run using the malignant outcome for 2011–2023 (Eurostat, 31 countries) and WHO data for both outcomes. It computed:

- deaths and age-standardised rates;
- a three-factor decomposition of the change between 2011–2013 and 2021–2023;
- constant-rate and trend projections, with a single-origin backtest;
- ratios of benign and uncertain to malignant CNS deaths by country and period.

The protocol was revised in response. The coding component and the composite outcome were introduced, trend extrapolation was subordinated to out-of-sample evaluation, years of life lost were removed as an outcome, and the UK was added. The primary analysis (section 5.1, Set L, 2003–2005 against 2020–2022, both outcomes) had not been run when this protocol was deposited.

## 8. Reporting, deviations and data sharing

Reporting follows GATHER (Guidelines for Accurate and Transparent Health Estimates Reporting) and RECORD. Deviations from this protocol will be listed with reasons in the publication. Code, filtered extracts and metadata will be archived with a DOI; raw WHO files are not redistributed and can be obtained from the WHO website.

## 9. Planned outputs

- Figure: recorded deaths and age-standardised rates, malignant and composite outcomes, Europe and regions.
- Figure: decomposition by country into demographic, coding and residual rate components.
- Figure: chained annual components with coding changes marked.
- Figure: projections to 2040.
- Table: deaths by country, both outcomes, baseline and end period.
- Table: out-of-sample performance of the projection models.
- Supplement: data audit, country sets, coding changes, sensitivity analyses, GATHER and RECORD checklists.

## References

Das Gupta P. Standardization and Decomposition of Rates: A User's Manual. Current Population Reports P23-186. Washington, DC: US Bureau of the Census; 1993.

Kuang Z, Tu J, Jiang M, et al. Global disease burden, trends, and inequalities of brain and central nervous system cancers, 1990–2021: a population-based study with projections to 2036. World Neurosurg 2025;198:123970. https://doi.org/10.1016/j.wneu.2025.123970

Møller B, Fekjær H, Hakulinen T, et al. Prediction of cancer incidence in the Nordic countries: empirical comparison of different approaches. Stat Med 2003;22:2751–2766. https://doi.org/10.1002/sim.1481

Riebler A, Held L. Projecting the future burden of cancer: Bayesian age-period-cohort analysis with integrated nested Laplace approximations. Biom J 2017;59:531–549. https://doi.org/10.1002/bimj.201500263

Yuan G, Miao D, Xu B, et al. Global, regional, and national burden of brain and other central nervous system cancers from 1990 to 2021, and projections to 2035: a systematic analysis for the global burden of disease study 2021. Int J Surg 2026;112:286–302. https://doi.org/10.1097/JS9.0000000000003497

## Appendix A. Documented changes in automated coding, 2011 onwards

Source: Eurostat country metadata on causes of death (table 18.3.2, 'Automated coding', in the metadata pages or their annexes), retrieved 30 October 2026, extracted with the project code (step 6) and checked against the source tables. A change is recorded when the automated-coding status or the family of the coding system differs from the previous documented year; software version updates are not counted. Belgium reports two regions separately. Luxembourg changed from STYX to Iris on 1 August 2011, so 2012 is the first full year with the new system. In total, 18 changes in 16 countries; no change between 2011 and 2023 was documented for the other countries.

| Country | First data year with change | Previous system | New system |
|---|---|---|---|
| Austria | 2015 | no automated coding | Iris |
| Belgium | 2012 | Flanders and Brussels: no automated coding | Flanders and Brussels: Iris |
| Belgium | 2017 | Wallonia: STYX | Wallonia: Iris |
| Cyprus | 2012 | MICAR | Iris |
| Denmark | 2015 | ACME | Iris, ACME |
| Finland | 2023 | own application | Iris |
| France | 2012 | Iris, ACME, MICAR, STYX, manual | Iris, ACME, MICAR, manual |
| France | 2016 | Iris, ACME, MICAR, manual | Iris, MUSE, manual |
| Iceland | 2012 | ACME | Iris |
| Ireland | 2017 | ACME, MICAR | Iris |
| Italy | 2016 | ACME, MICAR | Iris |
| Latvia | 2014 | ACME | Iris, ACME |
| Luxembourg | 2012 | Iris, STYX | Iris |
| Malta | 2013 | no automated coding | Iris |
| Netherlands | 2013 | no automated coding | Iris, manual |
| Spain | 2014 | no automated coding | Iris |
| Sweden | 2012 | ACME, MIKADO | Iris |
| Switzerland | 2023 | Iris | manual |

## Appendix B. Country sets

Derived from data availability only (WHO and Eurostat files listed in section 4). Mean deaths are C70–C72 deaths per year, 2017–2019.

| Country | Set L | Set L, complete case | Set W | Mean deaths 2017–2019 | Country-level results |
|---|---|---|---|---|---|
| Austria | yes | yes | yes | 567 | yes |
| Belgium | yes | yes | yes | 667 | yes |
| Bulgaria | no | no | yes | 667 | yes |
| Croatia | no | no | yes | 415 | yes |
| Cyprus | no | no | yes | 47 | no |
| Czechia | yes | yes | yes | 756 | yes |
| Denmark | yes | yes | yes | 479 | yes |
| Estonia | yes | yes | yes | 108 | yes |
| Finland | yes | yes | yes | 373 | yes |
| France | yes | yes | yes | 4,005 | yes |
| Germany | yes | no | yes | 6,075 | yes |
| Greece | no | no | yes | 1,063 | yes |
| Hungary | yes | yes | yes | 658 | yes |
| Iceland | yes | yes | yes | 23 | no |
| Ireland | no | no | yes | 295 | yes |
| Italy | yes | yes | yes | 4,292 | yes |
| Latvia | yes | yes | yes | 168 | yes |
| Lithuania | yes | yes | yes | 241 | yes |
| Luxembourg | yes | yes | yes | 31 | no |
| Malta | yes | yes | yes | 30 | no |
| Netherlands | yes | yes | yes | 1,046 | yes |
| Norway | yes | no | yes | 349 | yes |
| Poland | yes | yes | yes | 3,010 | yes |
| Portugal | no | no | yes | 893 | yes |
| Romania | yes | yes | yes | 1,609 | yes |
| Serbia | yes | yes | yes | 669 | yes |
| Slovakia | yes | yes | yes | 394 | yes |
| Slovenia | yes | yes | yes | 179 | yes |
| Spain | yes | yes | yes | 3,118 | yes |
| Sweden | yes | yes | yes | 654 | yes |
| Switzerland | yes | yes | yes | 565 | yes |
| United Kingdom | yes | yes | yes | 4,503 | yes |
