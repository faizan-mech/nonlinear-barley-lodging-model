#!/usr/bin/env bash
# One SLURM job of the V35 runs: bash JOB_V35.sh <stage> [nlanes]
# Production lanes take their lane number from SLURM_ARRAY_TASK_ID.
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Submit with sbatch (SUBMIT_V35.sh); do not run on the login node." >&2; exit 2; }
stage=${1:?stage required}
nlanes=${2:-1}
case "$stage" in calibration|conv_modal|conv_release|pilot|conv_records|conv_report|production|ensemble|postures|collect) ;; *) echo "bad stage $stage" >&2; exit 2;; esac
[[ "$nlanes" =~ ^[0-9]+$ ]] || exit 2
lane=${SLURM_ARRAY_TASK_ID:-1}
root=${V35_CODE_ROOT:?V35_CODE_ROOT not set; submit with SUBMIT_V35.sh}
[[ -f "$root/RUN_V35.m" ]] || { echo "Missing RUN_V35.m under $root" >&2; exit 2; }
export V35_OUT=${V35_OUT:?V35_OUT not set}
mkdir -p "$V35_OUT/Locks" "$V35_OUT/JobStorage"
exec 9>"$V35_OUT/Locks/${stage}_${lane}.lock"
flock -n 9 || { echo "Another job is writing $stage lane $lane; stopping." >&2; exit 73; }
export V35_JOBSTORAGE="$V35_OUT/JobStorage/${SLURM_JOB_ID}_${lane}"
mkdir -p "$V35_JOBSTORAGE"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
if ! type module >/dev/null 2>&1 && [[ -r /etc/profile.d/modules.sh ]]; then source /etc/profile.d/modules.sh; fi
if type module >/dev/null 2>&1; then module load ${V35_MATLAB_MODULE:-matlab/2022b}; fi
command -v matlab >/dev/null || { echo 'MATLAB is unavailable.' >&2; exit 127; }
cd -- "$root"
printf 'V35 root %s | out %s | stage %s | lane %s/%s | cpus %s\n' "$root" "$V35_OUT" "$stage" "$lane" "$nlanes" "${SLURM_CPUS_PER_TASK:-?}"
matlab -singleCompThread -batch "r=getenv('V35_CODE_ROOT'); cd(r); addpath(r,'-begin'); RUN_V35('$stage',$lane,$nlanes)"
rm -rf "$V35_JOBSTORAGE"
