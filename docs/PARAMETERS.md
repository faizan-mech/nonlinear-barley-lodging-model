# Reference parameters

Physical inputs are defined in `model/core/v35_params.m`; cases override them in `v35_cases.m` and `v35_apply_variant.m`.

| Quantity | Value or construction |
| --- | --- |
| Culm length | 0.83 m |
| Outer radius below the taper | 1.8857 mm; diameter 3.7714 mm |
| Taper | Starts at 0.5L; outer radius decreases linearly to 0.57 of its untapered value |
| Middle and upper wall thickness | 0.390 mm |
| Lower wall thickness | 0.390 mm multiplied by a separate construction profile with anchors 0.075, 0.075 and 0.042 N m at s = 0, 0.20 and 0.315 m, divided by 0.042 N m |
| Basal wall thickness | Approximately 0.69643 mm from that construction |
| Section geometry | A = pi(Ro^2 - Ri^2), I = pi(Ro^4 - Ri^4)/4, Z = I/Ro |
| Young's modulus | Derived so the harmonic mean of EI on s = 0.315--0.515 m is 0.0246 N m^2; 4.515 GPa for the reference geometry |
| Poisson ratio | 0.30 |
| Shear factor | 0.55, representative of Cowper's (1966) values for the hollow culm sections |
| Culm-plus-leaf mass parameter | 3.110 g |
| Head mass | 1.875 g; total of the two mass parameters is 4.985 g |
| Head centre-of-mass offset | 41.2 mm |
| Head centre-of-pressure offset | 63.2 mm |
| Head central rotary inertia | 1.1342711191582338e-6 kg m^2 |
| Upright and end-on projected areas | 1.605e-3 and 9.39e-5 m^2 |
| Head drag coefficient | 0.866 |
| Culm drag coefficient / attenuation | 1.0 / exponential attenuation coefficient 2.6 |
| Capacity anchors | s = [0, 0.20, 0.315, 0.515, 0.63, 0.83] m; capacity = [0.075, 0.075, 0.042, 0.042, 0.021, 0.021] N m; linear interpolation |
| Support targets | 0.97 Hz free-standing; 1.36 Hz supported |
| Canopy window | Lower/upper centres 0.30L/0.625L; width 0.03L; cut-off 0.75L |
| Structural damping | 0.08 in the first ten gravity-equilibrium modes |
| Reference inflow | Iu = 0.40, Lu/L = 1, von Karman spectrum |
| Wind synthesis / ramp / record | 0.625 ms / 2 s / 60 s |
| Production mesh / nominal time step | 61 nodes / 1.25 ms |
| Root capacity law | 0.58 tau_s d^3; divide by 4.4 shoots when root sharing is enabled; d = 44 mm |

The lower-wall construction and the failure-capacity profile are separate: a capacity perturbation does not change the wall or the stiffness. Geometry perturbations keep the modulus derived from the reference geometry.

Support stiffnesses are calculated by `v35_calibrate`; `examples/reference_demo.m` exports them with the other model properties.
