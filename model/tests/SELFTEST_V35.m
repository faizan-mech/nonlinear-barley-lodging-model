function ok = SELFTEST_V35()
%SELFTEST_V35 Element, tangent, calibration and short-trajectory checks.
% Element regression against the reference element, finite-difference
% tangents of rigid-head inertia, head and culm drag, and supports,
% reference calibration, and one short record with failure localisation.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'core'), fullfile(root, 'tests'), fullfile(root, 'tests', 'ref'));
ok = true;
e = t_element_regression(true);           ok = ok && all(e < 1e-10);
e = t_tangent('crossflow', 'rotational'); ok = ok && e(1) < 1e-6 && e(4) < 1e-6;
e = t_tangent('projection', 'rotational');ok = ok && e(1) < 1e-6;
e = t_tangent('upright', 'winkler');      ok = ok && e(1) < 1e-6;
P = v35_params();
C = v35_calibrate(P, 'both', [], 241);
fprintf('calibration: k_theta %.4f N m/rad, k_M %.4f N/rad, f1 free %.4f, canopy %.4f Hz\n', ...
    C.k_theta, C.k_M, C.f_free(1), C.f_canopy(1));
ok = ok && abs(C.f_free(1) - P.f_free) < 1e-6 && abs(C.f_canopy(1) - P.f_canopy) < 1e-6;
P.k_theta = C.k_theta; P.k_M = C.k_M;
m = v35_prepare(v35_build(P, 31));
r = v35_record(m, P, 14, 1, struct('dt', 5e-3, 'T', 6));
fprintf('record U = 14 m/s, seed 1: %s, failed %d at t = %.3f s, s/L = %.3f, lambda = %.4f, budget closure %.4f\n', ...
    r.status, r.failed, r.t_fail, r.s_fail/m.L, r.lam_fail, r.budget.closure);
ok = ok && r.ok && r.failed && abs(r.lam_fail - 1) <= 2e-3 && abs(r.budget.closure - 1) < 5e-3;   % nodal lumping, O(h)
if ok, fprintf('SELFTEST_V35 passed\n'); else, fprintf('SELFTEST_V35 FAILED\n'); end
end
