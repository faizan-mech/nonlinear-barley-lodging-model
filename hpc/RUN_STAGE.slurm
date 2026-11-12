#!/usr/bin/env bash
# Portable launcher; submit from the repository root with sbatch.
# Usage: sbatch [resources] hpc/RUN_STAGE.slurm STAGE [LANES]
#SBATCH --job-name=barley
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%A_%a.out
set -euo pipefail

[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Submit this file with sbatch.' >&2; exit 2; }
stage=${1:?Expected a stage}
lanes=${2:-1}
lane=${SLURM_ARRAY_TASK_ID:-1}
[[ "$lanes" =~ ^[1-9][0-9]*$ && "$lane" =~ ^[1-9][0-9]*$ ]] || exit 2
((lane <= lanes)) || { echo 'Array task exceeds the declared lane count.' >&2; exit 2; }

root=${BARLEY_ROOT:-${SLURM_SUBMIT_DIR:?Submit from the repository root}}
root=$(cd -- "$root" && pwd -P)
[[ -f "$root/model/RUN_V35.m" ]] || { echo 'Cannot locate model/ under the repository root.' >&2; exit 2; }
export V35_OUT=${V35_OUT:-"$root/results/v35"}
[[ "$V35_OUT" = /* ]] || { echo 'Use an absolute output path.' >&2; exit 2; }

case "$stage" in
  calibration|conv_modal|conv_release|pilot|conv_records|conv_report|production|ensemble|postures|collect) ;;
  *) echo 'Unknown stage.' >&2; exit 2 ;;
esac

mkdir -p "$V35_OUT/Locks"
exec 9>"$V35_OUT/Locks/${stage}_${lane}.lock"
flock -n 9 || { echo 'Another job owns this stage/lane.' >&2; exit 73; }

temp_parent=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}
job_storage=$(mktemp -d "$temp_parent/barley_${SLURM_JOB_ID}_${lane}.XXXXXX")
trap 'rm -rf -- "$job_storage"' EXIT
export V35_JOBSTORAGE="$job_storage"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
if [[ -n "${MATLAB_MODULE:-}" ]]; then
  if ! type module >/dev/null 2>&1 && [[ -r /etc/profile.d/modules.sh ]]; then source /etc/profile.d/modules.sh; fi
  type module >/dev/null 2>&1 || { echo 'Module command unavailable; put MATLAB on PATH or initialize modules.' >&2; exit 127; }
  module load "$MATLAB_MODULE"
fi
command -v matlab >/dev/null || { echo 'MATLAB is unavailable.' >&2; exit 127; }

export BARLEY_STAGE="$stage" BARLEY_LANE="$lane" BARLEY_LANES="$lanes"
cd -- "$root/model"
printf 'Barley %s lane %s/%s, CPUs %s, output %s\n' "$stage" "$lane" "$lanes" "${SLURM_CPUS_PER_TASK:-?}" "$V35_OUT"
matlab -singleCompThread -batch "addpath(pwd,'-begin'); RUN_V35(getenv('BARLEY_STAGE'),str2double(getenv('BARLEY_LANE')),str2double(getenv('BARLEY_LANES')))"
