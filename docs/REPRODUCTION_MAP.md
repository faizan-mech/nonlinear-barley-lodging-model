# What this repository reproduces

| Paper component | Entry point / relevant code | Status |
| --- | --- | --- |
| Rod, head and loading equations | `model/core/v35_residual.m`, `v35_mass.m`, `v35_build.m`, `v35_sections.m` | Included |
| Reference inputs and sensitivity cases | `v35_params.m`, `v35_cases.m` | Included (66 cases) |
| Support calibration and modal structure | `RUN_V35('calibration')` | Included |
| Spatial and temporal convergence | `calibration`, `conv_modal`, `conv_release`, `pilot`, `conv_records`, `conv_report` stages | Included; see the order below |
| Per-seed failure thresholds | `RUN_V35('production', lane, nlanes)` | Included; 24 seeds by default |
| Fixed-speed ensembles and root/stem first-failure probabilities | `RUN_V35('ensemble')` | Included; requires completed threshold summaries |
| Reference sway trajectories continued after crossing | `RUN_V35('postures')` | Included; requires the reference threshold summary |
| Moment budgets | `v35_budget.m`, `v35_record.m` | Included |
| Small-deflection, steady-flow and oscillatory-flow benchmarks | Benchmark drivers | Not included |
| Figure generation | Figure scripts | Not included |
| Production result records | Per-record outputs | Not included; regenerate with the stages above |

Run the convergence stages in this order:

1. `calibration`.
2. `conv_modal`, `conv_release` and `pilot` (independent after calibration).
3. `conv_records` after the pilot.
4. `conv_report` after all convergence branches finish.

The production mesh and step are set in `V35_SETTINGS.m`; inspect the convergence report before changing them. Use a new output directory after changing any input.

