#!/usr/bin/env bash
# Wait for the main Mantel/MRM job, then run sensitivity, figures, validation,
# and the summary. Use only after 02_mantel_mrm.R has been started.
set -euo pipefail
cd /d/seed_meta_new/Downstream_analysis
RSCRIPT="/d/Program Files/R/R-4.5.2/bin/Rscript.exe"
LOG=result/ASV_host_phylo/logs
MAIN_LOG="$LOG/02_mantel_mrm.log"

while true; do
  if grep -q "DONE 02_mantel_mrm" "$MAIN_LOG"; then break; fi
  if grep -q "Error" "$MAIN_LOG"; then
    echo "02 failed; see $MAIN_LOG" >&2
    exit 1
  fi
  sleep 20
done

"$RSCRIPT" --vanilla result/ASV_host_phylo/scripts/03_sensitivity.R > "$LOG/03_sensitivity.log" 2>&1
"$RSCRIPT" --vanilla result/ASV_host_phylo/scripts/04_figures.R > "$LOG/04_figures.log" 2>&1
"$RSCRIPT" --vanilla result/ASV_host_phylo/scripts/05_validate.R > "$LOG/05_validate.log" 2>&1
"$RSCRIPT" --vanilla result/ASV_host_phylo/scripts/06_summarize.R > "$LOG/06_summarize.log" 2>&1
echo "REMAINING PIPELINE DONE"
