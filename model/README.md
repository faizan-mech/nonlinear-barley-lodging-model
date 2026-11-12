# model/

MATLAB source of the V35 barley shoot model.

- `RUN_V35.m`, `V35_SETTINGS.m`: stage driver and run settings of the model.
- `core/`: rod, head, supports, wind synthesis, time integration, failure criterion, threshold search and sensitivity cases. Physical inputs are in `core/v35_params.m`; cases are in `core/v35_cases.m`.
- `tests/`: self-test (`SELFTEST_V35`) and tangent and element checks.
- `hpc/`: SLURM submission scripts with automatic stage dependencies.

Instructions are in the repository [README](../README.md) and in [docs/](../docs/).
