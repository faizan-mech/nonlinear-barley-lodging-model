function R = reference_demo(outDir)
%REFERENCE_DEMO Serial reference calibration and a short trajectory.
% Run from the repository root: addpath('examples'); R = reference_demo;
% Integrates a 3 s trajectory at 5 m/s.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'model'), fullfile(root, 'model', 'core'), '-begin');
if nargin < 1, outDir = fullfile(root, 'results', 'demo'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end

S = V35_SETTINGS();
% Fresh calibration cache.
cacheDir = tempname; mkdir(cacheDir);
cleanupCache = onCleanup(@()rmdir(cacheDir, 's')); %#ok<NASGU>
[m, P, C] = v35_case_model('ref', S.N, cacheDir, S.Ncal);

fprintf('Reference: %d nodes, dt %.6g s, shear factor %.3g\n', ...
    m.n, S.dt, P.k_shear);
fprintf('E %.9g Pa; k_theta %.9g N m/rad; k_M %.9g N/rad\n', ...
    m.E, C.k_theta, C.k_M);
fprintf('Free/supported calibration frequencies: %.6f / %.6f Hz\n', ...
    C.f_free(1), C.f_canopy(1));

% Synthesise the full 60 s record and integrate only its first 3 s.
opts = struct('dt', S.dt, 'T', 3, 'stop', true, 'budget', false, ...
    'hist_every', 8, 'tail_T', 0);
r = v35_record(m, P, 5.0, 1, opts);
assert(r.ok, 'Reference demonstration failed: %s', r.status);

R = struct('parameters', P, 'supports', C, 'E', m.E, ...
    'frequency_Hz', m.f_grav, 'record', r, 'demo_exposure_s', opts.T);
save(fullfile(outDir, 'reference_demo.mat'), 'R', '-v7');

names = {'L_m'; 'R_mid_m'; 't_mid_m'; 'E_Pa'; 'shear_factor'; ...
    'k_theta_Nm_per_rad'; 'k_M_N_per_rad'; 'nodes'; 'dt_s'; 'demo_exposure_s'};
values = [P.L; P.R_mid; P.t_mid; m.E; P.k_shear; C.k_theta; C.k_M; ...
    m.n; S.dt; opts.T];
writetable(table(names, values, 'VariableNames', {'parameter', 'value'}), ...
    fullfile(outDir, 'parameters.csv'));

H = double(r.hist);
assert(~isempty(H), 'No time history returned.');
T = array2table(H, 'VariableNames', {'time_s', 'tip_x_m', 'tip_y_m', ...
    'tip_theta_rad', 'stem_index', 'peak_node', 'base_moment_Nm', 'wind_x_ms'});
writetable(T, fullfile(outDir, 'trajectory.csv'));
fig = figure('Visible', 'off', 'Color', 'w');
cleanupFigure = onCleanup(@()close(fig)); %#ok<NASGU>
subplot(2,1,1); plot(H(:,1), H(:,2)/P.L, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Tip x/L'); grid on;
title('Short demonstration: U = 5 m/s, seed 1');
subplot(2,1,2); plot(H(:,1), H(:,5), 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Maximum stem failure index'); grid on;
print(fig, fullfile(outDir, 'reference_demo.png'), '-dpng', '-r150');
fprintf('Demo completed; outputs: %s\n', outDir);
end
