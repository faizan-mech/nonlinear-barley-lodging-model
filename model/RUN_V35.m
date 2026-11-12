function RUN_V35(stage, lane, nlanes)
%RUN_V35 Stage driver of the V35 barley shoot model.
%
%   RUN_V35('calibration')        supports of every case, frequency tables,
%                                  modal evolution with wind
%   RUN_V35('conv_modal')         frequencies and mean-wind equilibrium vs mesh
%   RUN_V35('conv_release')       free release vs step and mesh
%   RUN_V35('pilot')              reference U50 on the pilot discretisation
%   RUN_V35('conv_records')       30 paired records at the pilot U50, all cells
%   RUN_V35('conv_report')        acceptance table and recommended N, dt
%   RUN_V35('production', k, K)   threshold cases of lane k of K (all cases, or
%                                  only the ids in the environment variable
%                                  V35_CASES, comma-separated)
%   RUN_V35('postures')           swept centreline region of the 24 reference
%                                  records at U50
%   RUN_V35('ensemble')           sway ensembles and root-stem competition grid
%   RUN_V35('collect')            U50 table of all cases (CSV)
%
% Set the production N and dt in V35_SETTINGS after 'conv_report'.

if nargin < 2, lane = 1; end
if nargin < 3, nlanes = 1; end
if ischar(lane), lane = str2double(lane); end
if ischar(nlanes), nlanes = str2double(nlanes); end
S = V35_SETTINGS();
addpath(fullfile(S.root, 'core'));
if ~exist(S.outDir, 'dir'), mkdir(S.outDir); end
fprintf('RUN_V35 %s (lane %d of %d) on %s\n', stage, lane, nlanes, datestr(now));
switch stage
    case 'calibration',  stage_calibration(S);
    case 'conv_modal',   stage_conv_modal(S);
    case 'conv_release', start_pool(); stage_conv_release(S);
    case 'pilot',        start_pool(); stage_pilot(S);
    case 'conv_records', start_pool(); stage_conv_records(S);
    case 'conv_report',  stage_conv_report(S);
    case 'production',   start_pool(); stage_production(S, lane, nlanes);
    case 'ensemble',     start_pool(); stage_ensemble(S);
    case 'postures',     start_pool(); stage_postures(S);
    case 'collect',      stage_collect(S);
    otherwise, error('RUN_V35:stage', 'Unknown stage %s.', stage);
end
fprintf('RUN_V35 %s finished on %s\n', stage, datestr(now));
end

%% ======================================================================
function stage_calibration(S)
d = fullfile(S.outDir, 'calibration'); if ~exist(d, 'dir'), mkdir(d); end
C = v35_cases();
ids = [{'ref'}, unique({C(~strcmp({C.calib}, 'none')).id})];
T = struct('sid', {}, 'S', {});
for i = 1:numel(ids)
    cs = v35_find_case(ids{i});
    P = v35_apply_variant(v35_params(), cs.ov);
    Sx = v35_supports(cs, P, S.cacheDir, S.Ncal);
    T(end+1) = struct('sid', ids{i}, 'S', Sx); %#ok<AGROW>
end
fid = fopen(fullfile(d, 'calibration.txt'), 'w');
R = T(1).S;
out = @(varargin) fprintf_both(fid, varargin{:});
out('V35 calibration (mesh %d nodes)\n', S.Ncal);
out('Reference: clamped f1-f3 = %.4f %.4f %.4f Hz\n', R.f_clamp);
out('  k_theta = %.4f N m/rad for f1 = 0.97 Hz (0.945-0.995 Hz: %.4f-%.4f); Crook & Ennos %.2f-%.2f: inside = %d\n', ...
    R.k_theta, R.k_theta_interval, R.crook_ennos, R.in_crook_ennos);
out('  halving k_theta changes f1 by %.2f %%; a 20 %% drop needs k_theta x %.3f\n', 100*R.halving_drop, R.k_theta_ratio_for_20pc_drop);
out('  free-standing f1-f3 = %.4f %.4f %.4f Hz\n', R.f_free);
out('  k_M = %.4f N/rad for 1.36 Hz; in-canopy f1-f3 = %.4f %.4f %.4f Hz\n', R.k_M, R.f_canopy);
out('  k_M against the in-canopy target frequency:\n');
out('    f = %.2f Hz: k_M = %.4f\n', R.canopy_table.');
out('\nCase calibrations:\n%-12s %-5s %9s %9s %9s %8s %8s %8s\n', 'case', 'rule', 'k_theta', 'k_M', 'k_W', 'f1free', 'f1can', 'f2can');
for i = 1:numel(T)
    x = T(i).S;
    out('%-12s %-5s %9.4f %9.4f %9.4f %8.4f %8.4f %8.4f\n', T(i).sid, x.rule, x.k_theta, x.k_M, x.k_W, x.f_free(1), x.f_canopy(1), x.f_canopy(2));
end
% modal evolution with the mean wind (reference, 241 nodes)
[m, P] = v35_case_model('ref', 241, S.cacheDir, S.Ncal);
Pf = P; Pf.k_M = 0; mf = v35_prepare(v35_build(Pf, 241));
[f0, Phi0] = v35_modes(m, m.q_grav, [0 0], 3);
[ff, Phif] = v35_modes(mf, mf.q_grav, [0 0], 3);
nU = numel(S.modal_U); E = struct('U', S.modal_U, 'Ueq', nan(1,nU), 'f', nan(3,nU), 'tip', nan(2,nU), ...
    'lam', nan(1,nU), 's_lam', nan(1,nU), 'Mbase', nan(1,nU), 'phi1', nan(m.nf, nU), 'q', nan(m.nd, nU));
q = m.q_grav;
for i = 1:nU
    Ueq = v35_ueq(S.modal_U(i), P.wind.Iu*S.modal_U(i));
    [q, inf0] = v35_equilibrium(m, Ueq, q);
    if ~inf0.ok, break; end
    [f, Phi] = v35_modes(m, q, [Ueq 0], 3);
    [Mn, Mb] = v35_moments(m, q); [l, il] = max(abs(Mn)./m.BS);
    E.Ueq(i) = Ueq; E.f(:,i) = f; E.tip(:,i) = q(3*m.n-[2 1]); E.lam(i) = l; E.s_lam(i) = m.s(il);
    E.Mbase(i) = Mb; E.phi1(:,i) = Phi(:,1); E.q(:,i) = q;
end
out('\nModal evolution with the mean wind (reference, Iu = %.2f):\n   U     Ueq    f1     f2     f3    tip x/L  lam_mean  s/L\n', P.wind.Iu);
out('%5.1f %6.2f %6.3f %6.3f %6.3f %7.3f %8.3f %6.3f\n', [E.U; E.Ueq; E.f; E.tip(1,:)/m.L; E.lam; E.s_lam/m.L]);
fclose(fid);
modes = struct('s', m.s, 'canopy_f', f0, 'canopy_Phi', Phi0, 'free_f', ff, 'free_Phi', Phif, 'free', m.free);
v35_save_atomic(fullfile(d, 'calibration.mat'), struct('T', T, 'E', E, 'modes', modes));
end

function fprintf_both(fid, varargin)
fprintf(varargin{:}); fprintf(fid, varargin{:});
end

%% ======================================================================
function stage_conv_modal(S)
d = fullfile(S.outDir, 'convergence'); if ~exist(d, 'dir'), mkdir(d); end
Ns = S.conv.N_modal; nN = numel(Ns);
R = struct('N', Ns, 'f_canopy', nan(5,nN), 'f_free', nan(5,nN), 'f_clamp', nan(3,nN), ...
    'eq', []);
nU = numel(S.conv.U_eq_check);
EQ = struct('U', S.conv.U_eq_check, 'tip', nan(2,nU,nN), 'lam', nan(nU,nN), 's_lam', nan(nU,nN), ...
    'Mbase', nan(nU,nN), 'f1', nan(nU,nN), 'x31', nan(31,nU,nN));
s31 = v35_mesh(v35_params(), 31);
for k = 1:nN
    [m, P] = v35_case_model('ref', Ns(k), S.cacheDir, S.Ncal);
    R.f_canopy(:,k) = m.f_grav(1:5);
    Pf = P; Pf.k_M = 0; mf = v35_prepare(v35_build(Pf, Ns(k))); R.f_free(:,k) = mf.f_grav(1:5);
    Pc = P; Pc.clamped = true; Pc.k_M = 0; mc = v35_prepare(v35_build(Pc, Ns(k))); R.f_clamp(:,k) = mc.f_grav(1:3);
    q = m.q_grav;
    for i = 1:nU
        Ueq = v35_ueq(EQ.U(i), P.wind.Iu*EQ.U(i));
        q = v35_equilibrium(m, Ueq, q);
        [Mn, Mb] = v35_moments(m, q); [l, il] = max(abs(Mn)./m.BS);
        EQ.tip(:,i,k) = q(3*m.n-[2 1]); EQ.lam(i,k) = l; EQ.s_lam(i,k) = m.s(il); EQ.Mbase(i,k) = Mb;
        f = v35_modes(m, q, [Ueq 0], 1); EQ.f1(i,k) = f(1);
        [~, ia] = ismember(round(s31*1e9), round(m.s*1e9));   % nested meshes
        EQ.x31(:,i,k) = q(3*ia-2);
    end
    fprintf('N = %4d: canopy f1 %.6f, free f1 %.6f, clamp f1 %.6f Hz\n', Ns(k), R.f_canopy(1,k), R.f_free(1,k), R.f_clamp(1,k));
end
R.eq = EQ;
ref = nN;
R.err_f_canopy = R.f_canopy./R.f_canopy(:,ref) - 1;
R.err_f_free = R.f_free./R.f_free(:,ref) - 1;
R.err_lam = EQ.lam./EQ.lam(:,ref) - 1;
R.err_tip = squeeze(EQ.tip(1,:,:))./EQ.tip(1,:,ref).' - 1;
R.err_centreline = squeeze(max(abs(EQ.x31 - EQ.x31(:,:,ref)), [], 1))/v35_params().L;
R.err_f1_eq = EQ.f1./EQ.f1(:,ref) - 1;
v35_save_atomic(fullfile(d, 'conv_modal.mat'), struct('R', R));
fid = fopen(fullfile(d, 'conv_modal.txt'), 'w');
fprintf_both(fid, 'Relative error against N = %d\n  N    f1can     f2can     f3can     f1free    lam_eq(max over U)  tip_eq   centreline/L  f1_eq\n', Ns(ref));
for k = 1:nN
    fprintf_both(fid, '%4d %9.2e %9.2e %9.2e %9.2e %12.2e %12.2e %10.2e %10.2e\n', Ns(k), R.err_f_canopy(1:3,k), R.err_f_free(1,k), ...
        max(abs(R.err_lam(:,k))), max(abs(R.err_tip(:,k))), max(R.err_centreline(:,k)), max(abs(R.err_f1_eq(:,k))));
end
fclose(fid);
end

%% ======================================================================
function stage_conv_release(S)
d = fullfile(S.outDir, 'convergence'); if ~exist(d, 'dir'), mkdir(d); end
Ns = S.conv.N_release; dts = S.conv.dt_release;
[I, J] = ndgrid(1:numel(Ns), 1:numel(dts)); I = I(:); J = J(:);
nc = numel(I); res = cell(nc,1);
cacheDir = S.cacheDir; Ncal = S.Ncal;
parfor k = 1:nc
    res{k} = release_cell(Ns(I(k)), dts(J(k)), cacheDir, Ncal);
end
R = [res{:}];
iref = find([R.N] == max(Ns) & [R.dt] == min(dts));
fid = fopen(fullfile(d, 'conv_release.txt'), 'w');
fprintf_both(fid, 'Free release in still air, aerodynamic loads off (structural numerics only)\n');
fprintf_both(fid, 'small amplitude: frequency and damping against the damped eigenvalue; large amplitude: against N = %d, dt = %.4g ms\n', R(iref).N, 1e3*R(iref).dt);
fprintf_both(fid, '   N   dt(ms)  f_small err   zeta err    f_large err  peak-lam err  (aero on: f, zeta)\n');
for k = 1:nc
    x = R(k);
    fprintf_both(fid, '%4d %7.3f %11.2e %10.2e %12.2e %12.2e   (%.4f, %.4f)\n', x.N, 1e3*x.dt, x.small.f/x.fd - 1, x.small.zeta/x.zeta - 1, ...
        x.large.f/R(iref).large.f - 1, x.large_lam/R(iref).large_lam - 1, x.aero.f, x.aero.zeta);
end
fclose(fid);
v35_save_atomic(fullfile(d, 'conv_release.mat'), struct('R', R));
end

function x = release_cell(N, dt, cacheDir, Ncal)
[m, P] = v35_case_model('ref', N, cacheDir, Ncal);
f1 = m.f_grav(1); x = struct('N', N, 'dt', dt, 'f1', f1, 'zeta', P.zeta, 'fd', f1*sqrt(1 - P.zeta^2));
W = v35_wind_still(5);
o = struct('dt', dt, 'stop', false, 'hist_every', 1);
xg = m.q_grav(3*m.n-2);
% aerodynamics on: information only (still-air drag adds damping)
q0 = v35_equilibrium(m, 1.5);
out = v35_newmark(m, q0, W, o); x.aero = v35_decay(double(out.hist(:,1)), double(out.hist(:,2)), xg);
% aerodynamics off
ma = m; ma.cst = 0*m.cst; ma.ear.sa = 0;
out = v35_newmark(ma, q0, W, o); x.small = v35_decay(double(out.hist(:,1)), double(out.hist(:,2)), xg);
q1 = v35_equilibrium(m, 8);
out = v35_newmark(ma, q1, W, o); x.large = v35_decay(double(out.hist(:,1)), double(out.hist(:,2)), xg);
x.large_lam = out.lam_stem_max; x.large_tipmin = min(double(out.hist(:,2)));
end

%% ======================================================================
function stage_pilot(S)
d = fullfile(S.outDir, 'pilot');
set = struct('N', S.conv.pilot_N, 'dt', S.conv.pilot_dt, 'seeds', S.seeds, 'tol', S.tol, ...
    'outDir', d, 'cacheDir', S.cacheDir, 'Ncal', S.Ncal, 'max_hours', S.max_hours, ...
    'logfile', fullfile(d, 'pilot.log'));
if ~exist(d, 'dir'), mkdir(d); end
R = v35_threshold({'ref'}, set);
if R(1).complete
    pilot = R(1).stats; v35_save_atomic(fullfile(d, 'pilot.mat'), struct('pilot', pilot));
    fprintf('Pilot U50 = %.3f m/s [%.3f, %.3f]\n', pilot.U50, pilot.U50_ci);
end
end

%% ======================================================================
function stage_conv_records(S)
L = load(fullfile(S.outDir, 'pilot', 'pilot.mat')); U = L.pilot.U50;
d = fullfile(S.outDir, 'conv_records');
jobs = struct('cid', {}, 'N', {}, 'dt', {}, 'U', {}, 'seed', {}, 'file', {}, 'ropts', {});
for c = 1:size(S.conv.cells, 1)
    N = S.conv.cells(c,1); dt = S.conv.cells(c,2);
    for sd = S.conv.paired_seeds
        jobs(end+1) = struct('cid', 'ref', 'N', N, 'dt', dt, 'U', U, 'seed', sd, ...
            'file', fullfile(d, sprintf('N%03d_dt%06.3fms', N, 1e3*dt), sprintf('s%03d.mat', sd)), ...
            'ropts', struct('stop', false, 'budget', false)); %#ok<AGROW>
    end
end
set = struct('cacheDir', S.cacheDir, 'Ncal', S.Ncal, 'max_hours', S.max_hours, 'logfile', fullfile(d, 'conv_records.log'));
if ~exist(d, 'dir'), mkdir(d); end
v35_ensemble(jobs, set);
end

%% ======================================================================
function stage_conv_report(S)
d = fullfile(S.outDir, 'convergence');
Lm = load(fullfile(d, 'conv_modal.mat')); Rm = Lm.R;
Lr = load(fullfile(d, 'conv_release.mat')); Rr = Lr.R;
Lp = load(fullfile(S.outDir, 'pilot', 'pilot.mat')); U = Lp.pilot.U50;
cells = S.conv.cells; nc = size(cells,1);
dr = fullfile(S.outDir, 'conv_records');
seeds = S.conv.paired_seeds; ns = numel(seeds);
A = nan(nc, ns, 5);             % lam_max, s_max, first passage, t of max, lam_rms at site
for c = 1:nc
    for j = 1:ns
        f = fullfile(dr, sprintf('N%03d_dt%06.3fms', cells(c,1), 1e3*cells(c,2)), sprintf('s%03d.mat', seeds(j)));
        if exist(f, 'file')
            L = load(f); r = L.r;
            [~, is] = max(r.lam_env);
            A(c,j,:) = [r.lam_stem_max, r.s_stem_max, r.t_first_passage, r.t_stem_max, double(r.stats.lam_rms(is))];
        end
    end
end
[~, iref] = max(cells(:,1)./cells(:,2));
P = v35_params(); T = P.wind.T;
fid = fopen(fullfile(d, 'conv_report.txt'), 'w');
fprintf_both(fid, 'Paired records at the pilot U50 = %.3f m/s (%d seeds), reference cell N = %d, dt = %.4g ms\n', U, ns, cells(iref,1), 1e3*cells(iref,2));
fprintf_both(fid, '   N  dt(ms) |dlam|med  |dlam|p90  |dlam|max  |drms|med  outcome  site  f1err     rel_f     rel_zeta   pass\n');
best = []; bestCost = inf;
for c = 1:nc
    N = cells(c,1); dt = cells(c,2);
    dl = abs(A(c,:,1)./A(iref,:,1) - 1);
    fa = A(c,:,3) < T; fr = A(iref,:,3) < T;
    fa(isnan(A(c,:,3))) = false; fr(isnan(A(iref,:,3))) = false;
    outc = sum(fa == fr);
    hc = max(diff(v35_mesh(P, N)));   % largest element of the candidate mesh
    site = sum(abs(A(c,:,2) - A(iref,:,2)) <= hc + 1e-12);
    kN = find(Rm.N == N, 1); f1e = abs(Rm.err_f_canopy(1,kN));
    kr = find([Rr.N] == N & abs([Rr.dt] - dt) < 1e-12, 1);
    if isempty(kr)
        [~, kr] = min(abs([Rr.N] - N) + 1e3*abs([Rr.dt] - dt));
    end
    ef = abs(Rr(kr).small.f/Rr(kr).fd - 1); ez = abs(Rr(kr).small.zeta/Rr(kr).zeta - 1);
    p90 = sort(dl(~isnan(dl))); if isempty(p90), p90 = NaN; else, p90 = p90(max(1, ceil(0.9*numel(p90)))); end
    pass = f1e <= S.acc.f1 && ef <= S.acc.release_f && ez <= S.acc.release_zeta && ...
           p90 <= S.acc.lam && outc >= S.acc.outcome && site >= S.acc.site;
    drms = abs(A(c,:,5)./A(iref,:,5) - 1);
    fprintf_both(fid, '%4d %7.3f %10.2e %10.2e %10.2e %10.2e %6d/%d %4d %9.2e %9.2e %9.2e   %d\n', N, 1e3*dt, ...
        median(dl(~isnan(dl))), p90, max(dl), median(drms(~isnan(drms))), outc, ns, site, f1e, ef, ez, pass);
    cost = N/dt;
    if pass && c ~= iref && cost < bestCost, best = [N dt]; bestCost = cost; end
end
if isempty(best)
    fprintf_both(fid, '\nNo candidate meets every criterion; use the finest passing cell or refine.\n');
else
    fprintf_both(fid, '\nRecommended production discretisation: N = %d nodes, dt = %.4g ms (set in V35_SETTINGS)\n', best(1), 1e3*best(2));
end
fclose(fid);
v35_save_atomic(fullfile(d, 'conv_report.mat'), struct('A', A, 'cells', cells, 'best', best, 'U', U));
end

%% ======================================================================
function stage_production(S, lane, nlanes)
C = v35_cases(); ids = {C.id};
only = strtrim(strsplit(getenv('V35_CASES'), ','));
only = only(~cellfun(@isempty, only));
if ~isempty(only)
    bad = setdiff(only, ids);
    if ~isempty(bad), error('RUN_V35:cases', 'Unknown case ids: %s', strjoin(bad, ' ')); end
    ids = only;
end
mine = ids(mod((1:numel(ids)) - 1, nlanes) == lane - 1);
d = fullfile(S.outDir, 'production'); if ~exist(d, 'dir'), mkdir(d); end
U50 = [];
fr = fullfile(d, 'ref', 'summary.mat'); fp = fullfile(S.outDir, 'pilot', 'pilot.mat');
if exist(fr, 'file'), L = load(fr); U50 = L.st.U50;
elseif exist(fp, 'file'), L = load(fp); U50 = L.pilot.U50; end
set = struct('N', S.N, 'dt', S.dt, 'seeds', S.seeds, 'tol', S.tol, 'outDir', d, ...
    'cacheDir', S.cacheDir, 'Ncal', S.Ncal, 'U50_ref', U50, 'max_hours', S.max_hours, ...
    'ropts_ref', S.ropts_ref, 'logfile', fullfile(d, sprintf('lane%02d.log', lane)));
fprintf('Lane %d of %d: %s\n', lane, nlanes, strjoin(mine, ' '));
v35_threshold(mine, set);
end

%% ======================================================================
function stage_ensemble(S)
d = fullfile(S.outDir, 'production'); de = fullfile(S.outDir, 'ensemble');
jobs = struct('cid', {}, 'N', {}, 'dt', {}, 'U', {}, 'seed', {}, 'file', {}, 'ropts', {});
L = load(fullfile(d, 'ref', 'summary.mat')); U50 = L.st.U50;
for f = S.ens.sway_factors
    for sd = S.seeds
        jobs(end+1) = struct('cid', 'ref', 'N', S.N, 'dt', S.dt, 'U', round(f*U50*1e5)/1e5, 'seed', sd, ...
            'file', fullfile(de, 'sway_ref', sprintf('f%4.2f_s%03d.mat', f, sd)), ...
            'ropts', struct('stop', false, 'hist_every', 4, 'budget', false)); %#ok<AGROW>
    end
end
fr = fullfile(d, 'root_stem', 'summary.mat');
if exist(fr, 'file')
    L = load(fr); Ur = L.st.U50;
    for f = S.ens.root_grid
        for sd = S.seeds
            jobs(end+1) = struct('cid', 'root_stem', 'N', S.N, 'dt', S.dt, 'U', round(f*Ur*1e5)/1e5, 'seed', sd, ...
                'file', fullfile(de, 'root_stem_grid', sprintf('f%4.2f_s%03d.mat', f, sd)), ...
                'ropts', struct('stop', true)); %#ok<AGROW>
        end
    end
end
set = struct('cacheDir', S.cacheDir, 'Ncal', S.Ncal, 'max_hours', S.max_hours, 'logfile', fullfile(de, 'ensemble.log'));
if ~exist(de, 'dir'), mkdir(de); end
v35_ensemble(jobs, set);
end

%% ======================================================================
function stage_postures(S)
% Swept centreline region at the reference U50: the 24 seeds of the sway
% ensemble at 1.00 U50, integrated without stopping, snapshots every 20 ms.
d = fullfile(S.outDir, 'production'); de = fullfile(S.outDir, 'ensemble', 'postures_ref');
if ~exist(de, 'dir'), mkdir(de); end
L = load(fullfile(d, 'ref', 'summary.mat')); U = round(L.st.U50*1e5)/1e5;
[m, P] = v35_case_model('ref', S.N, S.cacheDir, S.Ncal);
seeds = S.seeds; dt = S.dt;
files = arrayfun(@(sd) fullfile(de, sprintf('s%03d.mat', sd)), seeds, 'UniformOutput', false);
todo = find(~cellfun(@(f) exist(f, 'file') == 2, files));
fprintf('Postures at U = %.5f m/s: %d of %d records to run\n', U, numel(todo), numel(seeds));
parfor q = 1:numel(todo)
    k = todo(q);
    r = v35_postures(m, P, U, seeds(k), struct('dt', dt, 'every', 0.02));
    v35_save_atomic(files{k}, struct('r', r));
end
fprintf('Postures written to %s\n', de);
end

%% ======================================================================
function stage_collect(S)
d = fullfile(S.outDir, 'production');
C = v35_cases();
fid = fopen(fullfile(S.outDir, 'U50_table.csv'), 'w');
fprintf(fid, 'id,group,description,U50,U50_lo,U50_hi,U50_lognormal,U10,U90,change_vs_ref_pct,median_s_fail_over_L,n_seeds,n_nonmonotone,n_excluded,f1_Hz,k_theta,k_M,k_W\n');
ref = NaN;
fr = fullfile(d, 'ref', 'summary.mat');
if exist(fr, 'file'), L = load(fr); ref = L.st.U50; end
for i = 1:numel(C)
    f = fullfile(d, C(i).id, 'summary.mat');
    if ~exist(f, 'file'), fprintf(fid, '%s,%s,"%s",,,,,,,,,,,,,,,\n', C(i).id, C(i).group, C(i).desc); continue; end
    L = load(f); st = L.st;
    fprintf(fid, '%s,%s,"%s",%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%.4f,%d,%d,%d,%.5f,%.5f,%.5f,%.5f\n', ...
        C(i).id, C(i).group, C(i).desc, st.U50, st.U50_ci, st.U50_lognormal, st.U10, st.U90, ...
        100*(st.U50/ref - 1), median(st.s_fail(~isnan(st.s_fail)))/st.L, st.n, st.n_nonmonotone, st.n_excluded, ...
        st.f_grav(1), st.k_theta, st.k_M, st.k_W);
end
fclose(fid);
fprintf('Wrote %s\n', fullfile(S.outDir, 'U50_table.csv'));
end

%% ======================================================================
function start_pool()
if exist('OCTAVE_VERSION', 'builtin') ~= 0, return; end
p = gcp('nocreate');
if ~isempty(p), return; end
n = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(n) || n < 1, n = feature('numcores'); end
try
    cl = parcluster('Processes');
catch
    cl = parcluster('local');
end
jd = getenv('V35_JOBSTORAGE');
if isempty(jd), jd = fullfile(tempdir, sprintf('v35_%s_%d', getenv('SLURM_JOB_ID'), feature('getpid'))); end
if ~exist(jd, 'dir'), mkdir(jd); end
cl.JobStorageLocation = jd;
cl.NumWorkers = n;
p = parpool(cl, n);
parfevalOnAll(p, @maxNumCompThreads, 0, 1);
fprintf('Pool: %d workers, job storage %s\n', p.NumWorkers, jd);
end
