# Nonlinear barley lodging model

MATLAB code accompanying:

> M. Faizan, N. Freidoonimehr, M. R. Tucker and M. Arjomandi, *Nonlinear modelling of barley shoot reconfiguration and local stem failure under turbulent wind*, submitted to the Journal of Fluids and Structures (2026).

A planar, geometrically exact rod model of a barley shoot with a rigid head, compliant root anchorage and effective canopy support. Turbulent wind drives the response, and local bending and anchorage capacities define stem and root failure. Head properties use the `ear` name in the code.

## Requirements

- MATLAB R2022b or later. Production runs require Parallel Computing Toolbox; the example runs serially.
- Bash, SLURM and `flock` for cluster runs.

## Quick start

From the repository root in MATLAB:

```matlab
addpath('model', 'model/core', 'model/tests', 'examples');
assert(SELFTEST_V35(), 'V35 self-test failed');
R = reference_demo;
```

The example calibrates the supports and runs a 3 s trajectory at 5 m/s. Parameters, the trajectory and a plot are saved in `results/demo/`.

## Failure thresholds

```matlab
addpath('model');
setenv('V35_OUT', fullfile(pwd, 'results', 'v35'));
RUN_V35('calibration');
setenv('V35_CASES', 'ref');
RUN_V35('production', 1, 1);
RUN_V35('collect');
```

Results are written to `results/v35/U50_table.csv`. Set `V35_CASES` to comma-separated case ids, or use `setenv('V35_CASES', '')` for all 66 cases defined in [v35_cases.m](model/core/v35_cases.m).

Defaults are 61 nodes, a 1.25 ms step, 24 seeds and 60 s records. The threshold search uses a 0.5 % bracket and reports observed non-monotonicity. Use a new output directory after changing inputs or discretisation.

## Documentation

- [Parameters](docs/PARAMETERS.md)
- [Analysis stages and included components](docs/REPRODUCTION_MAP.md)
- [SLURM runs](docs/HPC.md)

The release contains the model and analysis code. Benchmark drivers, manuscript figure scripts and production result files are not included.

## Citation and licence

See [CITATION.cff](CITATION.cff) and the [MIT licence](LICENSE).

Contact: Muhammad Faizan, Adelaide University — muhammad.faizan@adelaide.edu.au.
