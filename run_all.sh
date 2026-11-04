#!/usr/bin/env bash
# Full analysis pipeline, from the source files in data/raw/ to the tables in output/tables/ and the figures in
# output/figures/. With --fetch, the source files are downloaded first (internet access needed); without it,
# they must already be in data/raw/ (see data/README.md).
set -euo pipefail
export OMP_NUM_THREADS=1 PYTHONHASHSEED=0 TZ=UTC   # single-threaded BLAS/OpenMP
cd "$(dirname "$0")/scripts"
if [[ "${1:-}" == "--fetch" ]]; then
  python3 01_fetch_who_mdb.py
  python3 02_fetch_eurostat.py
  python3 02b_fetch_europop2025.py
  python3 02b_fetch_unwpp_metadata.py
  python3 02c_fetch_registries.py
else
  python3 01_fetch_who_mdb.py --no-download
fi
python3 03_harmonise.py
python3 06_coding_metadata.py
python3 07_country_sets.py
python3 08_extract_who_detail.py
python3 09_build_panel.py
python3 10_coding_changes.py
python3 11_decomposition.py
Rscript 12_trend_model.R
for fe in 2010 2013 2016; do FIT_END=$fe Rscript 12_trend_model.R; done
python3 13_rates_joinpoint.py
python3 14_population_projections.py
python3 15_projections.py prepare
Rscript 16_bapc.R
python3 15_projections.py evaluate
python3 17_secondary.py
python3 17b_coding_evidence.py
python3 17c_forecast_checks.py
python3 17d_country_comparisons.py
python3 17e_registry_england.py
python3 17f_registry_belgium.py
python3 18_figures.py
mkdir -p ../output
{ python3 --version; python3 -c "import pandas, numpy, scipy, statsmodels, matplotlib; print('pandas', pandas.__version__); print('numpy', numpy.__version__); print('scipy', scipy.__version__); print('statsmodels', statsmodels.__version__); print('matplotlib', matplotlib.__version__)"; pdftotext -v 2>&1 | head -1; R --version | head -1; Rscript -e 'for (p in c("glmmTMB","TMB","INLA","BAPC","MASS","Matrix")) cat(p, as.character(packageVersion(p)), "\n")'; } > ../output/software_versions.txt 2>&1
