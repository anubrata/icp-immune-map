#!/usr/bin/env bash
# Run the whole analysis (steps 1-4) and build the website data in site/data/.
#
#   ./run_all.sh                          every dataset in config/
#   ./run_all.sh config/GSE282182.yaml    just one dataset
#
# Then preview the site:  uv run python -m http.server --directory site 8000
set -euo pipefail
cd "$(dirname "$0")"

if [ "$#" -eq 0 ]; then
  rm -rf site/data            # full rebuild: drop data from datasets that no longer have a config
  set -- config/*.yaml
fi

for config in "$@"; do
  echo "================ $config ================"
  echo "--- Step 1: load data";                 uv run analysis/step1_load_data.py "$config"
  echo; echo "--- Step 2: differential expression"; uv run analysis/step2_differential_expression.py "$config"
  echo; echo "--- Step 3: immune annotation & cell types"; uv run analysis/step3_immune_annotation.py "$config"
  echo; echo "--- Step 4: export for the website";   uv run analysis/step4_export_site_data.py "$config"
  echo
done

echo "Done. Preview the site: uv run python -m http.server --directory site 8000  ->  http://localhost:8000"
