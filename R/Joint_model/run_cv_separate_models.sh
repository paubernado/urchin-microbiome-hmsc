#!/bin/bash
# =============================================================================
# Cross-validation, in-sample fit and results of the two separate models
# (PA and conditional abundance) already fitted on venus.
#
# Runs the tasks one after the other, with low priority (nice 19):
#   PA    : cv_random -> cv_site -> results
#   ABUND : insample -> cv_random -> cv_site -> results
# Folds are fitted one after the other: 8 cores in total (one per chain).
#
# USAGE (on venus, from ~/data/data/hsmc, where HMSC_STEP2_cv_results.R is):
#   screen -S hmsc_cv                  # open a screen session
#   bash run_cv_separate_models.sh     # PA first, then abundance (or: pa / abund)
#   Ctrl+A then D                      # detach (the jobs keep running)
#   screen -r hmsc_cv                  # come back later to check
#
# Progress: logs/*.log   (e.g. tail -f logs/pa_cv_site.log)
# =============================================================================

set -o pipefail

CORES=8                                  # cores per R process (= number of chains)
NICE=19          # lowest CPU priority: other users' jobs go first
PAR_FOLDS=no     # no = folds one after the other (8 cores); yes = 4 folds at once (32 cores)
CV_SAMPLES=125   # samples per chain in each CV fold model (8 chains x 125 = 1000 draws);
                 # same thin and transient as the full model. 0 = same as full model.

PA=results_phy_PA_SD_ranlev
MPA=$PA/HMSC_PA_with_phylogeny_logDepth_siteRandom.rds
AB=results_phy_ABU_SD_ranlev
MAB=$AB/HMSC_ABUN_phylogeny_logDepth_siteRandom.rds

# ---- conda environment with Hmsc ----
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hmsc_env || { echo "Could not activate hmsc_env"; exit 1; }
set -u   # after conda activate (its scripts use unset variables)

# One thread per R process. Without this, the linear algebra library (BLAS)
# used by each MCMC chain starts one thread per core of the machine, so 8
# chains end up using all 128 cores. The R worker processes inherit these.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# ---- checks ----
[ -f HMSC_STEP2_cv_results.R ] || { echo "HMSC_STEP2_cv_results.R not found here"; exit 1; }
[ -f "$MPA" ] || { echo "Not found: $MPA"; exit 1; }
[ -f "$MAB" ] || { echo "Not found: $MAB"; exit 1; }
mkdir -p logs
echo "Cores available: $(nproc) | free memory (GB): $(free -g | awk '/Mem:/ {print $7}')"

# The in-sample fit of the PA model already exists: copy it with the name STEP 2 expects
[ -f "$PA/MF_pa_insample.rds" ] || cp "$PA/MF_PA_with_phylogeny_insample.rds" "$PA/MF_pa_insample.rds"

# run_task <log name> <model> <outdir> <model file> <task>
run_task () {
  local name=$1 model=$2 outdir=$3 mfile=$4 task=$5
  echo "$(date '+%F %T') START $name"
  if nice -n "$NICE" Rscript HMSC_STEP2_cv_results.R --model "$model" --outdir "$outdir" --model-file "$mfile" \
       --task "$task" --n-parallel "$CORES" --cv-samples "$CV_SAMPLES" --parallel-folds "$PAR_FOLDS" > "logs/$name.log" 2>&1; then
    echo "$(date '+%F %T') OK    $name"
  else
    echo "$(date '+%F %T') ERROR $name  (see logs/$name.log)"; return 1
  fi
}

# Tasks run ONE AFTER THE OTHER (one CV at a time, 8 cores, low priority
# with nice so that other users' jobs go first):
#   PA    : cv_random -> cv_site -> results
#   ABUND : insample -> cv_random -> cv_site -> results

WHICH=${1:-all}     # bash run_cv_separate_models.sh pa | abund | all
case "$WHICH" in pa|abund|all) ;; *) echo "Use: pa, abund or all"; exit 1;; esac

if [ "$WHICH" != abund ]; then
  run_task pa_cv_random pa "$PA" "$MPA" cv_random &&
  run_task pa_cv_site   pa "$PA" "$MPA" cv_site   &&
  run_task pa_results   pa "$PA" "$MPA" results
fi

if [ "$WHICH" != pa ]; then
  run_task abund_insample  abund "$AB" "$MAB" insample  &&
  run_task abund_cv_random abund "$AB" "$MAB" cv_random &&
  run_task abund_cv_site   abund "$AB" "$MAB" cv_site   &&
  run_task abund_results   abund "$AB" "$MAB" results
fi

wait
echo "$(date '+%F %T') All tasks finished. Results in $PA/ and $AB/"
