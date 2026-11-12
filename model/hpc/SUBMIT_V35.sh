#!/usr/bin/env bash
# Submit the V35 runs from the HPC login node.
#   bash model/hpc/SUBMIT_V35.sh calib-conv      calibration -> convergence study -> report
#   bash model/hpc/SUBMIT_V35.sh production [K]  all threshold cases on K nodes (default 4)
#   bash model/hpc/SUBMIT_V35.sh ensemble        sway ensembles and root-stem grid (after production)
#   bash model/hpc/SUBMIT_V35.sh extra K ids     only the listed case ids (comma-separated) on K nodes
#   bash model/hpc/SUBMIT_V35.sh postures        swept centreline region at the reference U50
#   bash model/hpc/SUBMIT_V35.sh collect         U50 table
#   bash model/hpc/SUBMIT_V35.sh status
# Environment: V35_OUT (results root; default <repository>/model/results), V35_PARTITIONS (SLURM
# partitions, comma-separated; default: the cluster default), V35_NODE_CPUS (CPUs per node for
# the full-node stages; default 72), V35_MEM (MB, overrides the per-stage default of 4 GB per
# CPU + 8 GB; full-node stages use the whole node), V35_MATLAB_MODULE (default matlab/2022b).
set -euo pipefail
mode=${1:-status}
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
[[ -f "$root/RUN_V35.m" ]] || { echo "Incomplete package: $root" >&2; exit 2; }
export V35_CODE_ROOT="$root"
export V35_OUT=${V35_OUT:-$root/results}
partitions=${V35_PARTITIONS:-}
full=${V35_NODE_CPUS:-72}; [[ "$full" =~ ^[1-9][0-9]*$ ]] || { echo "bad V35_NODE_CPUS" >&2; exit 2; }
mkdir -p "$V35_OUT/Logs"
if [[ "$mode" == status ]]; then
  squeue -u "$USER" -o '%.12i %.24j %.10T %.12M %R' | grep -E 'JOBID|V35' || true
  for f in calibration/calibration.txt convergence/conv_report.txt U50_table.csv; do
    [[ -f "$V35_OUT/$f" ]] && { echo "== $f"; head -40 "$V35_OUT/$f"; }
  done
  ls "$V35_OUT"/production/*/summary.mat 2>/dev/null | wc -l | xargs echo "completed production cases:"
  exit 0
fi
command -v sbatch >/dev/null || { echo 'Run on the HPC login node.' >&2; exit 127; }
submit() {   # stage cpus wall dependency array nlanes
  local stage=$1 cpus=$2 wall=$3 dep=$4 array=$5 nl=${6:-1}
  local mem=${V35_MEM:-$(( cpus>=full ? 0 : 4000*cpus + 8000 ))}   # MB; 0 = whole node for full-node stages
  local -a a=(--parsable --job-name="V35_${stage}" --nodes=1 --ntasks=1
              --cpus-per-task="$cpus" --mem="$mem" --time="$wall" --chdir="$root" --export=ALL
              --output="$V35_OUT/Logs/%x_%A_%a.out" --error="$V35_OUT/Logs/%x_%A_%a.err")
  [[ -z "$partitions" ]] || a+=(--partition="$partitions")
  [[ -z "$dep" ]] || a+=(--dependency="$dep")
  [[ -z "$array" ]] || a+=(--array="$array")
  local r; r=$(sbatch "${a[@]}" "$root/hpc/JOB_V35.sh" "$stage" "$nl"); r=${r%%;*}
  printf '%s\t%s\t%s\n' "$r" "$stage" "$(date -u +%FT%TZ)" >> "$V35_OUT/Submitted_Jobs.tsv"
  echo "Submitted $stage: $r" >&2; printf '%s' "$r"
}
case "$mode" in
  calib-conv)
    c=$(submit calibration 4 06:00:00 '' '')
    mo=$(submit conv_modal 4 12:00:00 "afterok:$c" '')
    re=$(submit conv_release 16 1-00:00:00 "afterok:$c" '')
    pi=$(submit pilot "$full" 3-00:00:00 "afterok:$c" '')
    pr=$(submit conv_records "$full" 3-00:00:00 "afterok:$pi" '')
    submit conv_report 2 02:00:00 "afterok:$pr:$mo:$re" '' >/dev/null ;;
  production)
    K=${2:-4}; [[ "$K" =~ ^[1-9][0-9]*$ ]] || exit 2
    echo "Production with N and dt from V35_SETTINGS.m:"; grep -E "^S\.(N|dt) " "$root/V35_SETTINGS.m"
    submit production "$full" 3-00:00:00 '' "1-$K" "$K" >/dev/null ;;
  extra)
    K=${2:-1}; [[ "$K" =~ ^[1-9][0-9]*$ ]] || exit 2
    export V35_CASES=${3:?usage: SUBMIT_V35.sh extra K id1,id2,...}
    echo "Extra cases on $K node(s): $V35_CASES"; grep -E "^S\.(N|dt) " "$root/V35_SETTINGS.m"
    submit production "$full" 3-00:00:00 '' "1-$K" "$K" >/dev/null ;;
  postures) submit postures 24 12:00:00 '' '' >/dev/null ;;
  ensemble) submit ensemble "$full" 3-00:00:00 '' '' >/dev/null ;;
  collect)  submit collect 2 01:00:00 '' '' >/dev/null ;;
  *) echo "Unknown mode $mode" >&2; exit 2 ;;
esac
