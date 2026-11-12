# SLURM runs

Submit from the repository root. Set the output path and the local MATLAB module:

```bash
export V35_OUT="$PWD/results/v35"
export MATLAB_MODULE=matlab/2022b
```

Adjust the path, module name, partition and account for the cluster. Leave `MATLAB_MODULE` unset if MATLAB is already on `PATH`. The examples use 72-CPU nodes; `--mem=0` requests all node memory.

## Production

```bash
sbatch --cpus-per-task=4 --mem=24G --time=06:00:00 hpc/RUN_STAGE.slurm calibration
```

After calibration succeeds:

```bash
sbatch --cpus-per-task=72 --mem=0 --time=3-00:00:00 --array=1-4 hpc/RUN_STAGE.slurm production 4
```

Each array element runs one lane on one node. Set `V35_CASES=ref,ear_upright` to select cases. Preserve case selection and lane layout on restart. Once production finishes:

```bash
sbatch --cpus-per-task=2 --mem=16G --time=01:00:00 hpc/RUN_STAGE.slurm collect
```

`ensemble` and `postures` also require completed production summaries. Use production resources for `ensemble`; `postures` uses 24 CPUs and 12 h in the automatic submitter.

For convergence, follow the [stage order](REPRODUCTION_MAP.md). The automatic submitter `bash model/hpc/SUBMIT_V35.sh calib-conv` sets these dependencies; other modes are listed in its header.

The drivers use one MATLAB worker per allocated CPU. The automatic submitter reads `V35_OUT`, `V35_PARTITIONS`, `V35_NODE_CPUS`, `V35_MEM` and `V35_MATLAB_MODULE`. Keep the source unchanged while a study is running.
