function S = V35_SETTINGS()
%V35_SETTINGS Run settings. Edit production mesh and step after convergence.
S = struct();
S.root      = fileparts(mfilename('fullpath'));
S.outDir    = getenv('V35_OUT');               % results (set by the job scripts)
if isempty(S.outDir), S.outDir = fullfile(S.root, 'results'); end
S.cacheDir  = fullfile(S.outDir, 'cache');     % calibration cache, shared by all jobs
S.Ncal      = 481;                             % calibration mesh (nodes)
S.seeds     = 1:24;                            % wind seeds of every threshold case
S.tol       = 0.005;                           % relative bracket of each per-seed threshold

% ---- production discretisation (selected by the convergence study) ----
S.N         = 61;                              % nodes (31, 61, 121, 241)
S.dt        = 1.25e-3;                         % s

% ---- convergence study ----
S.conv.N_modal   = [31 61 121 241 481];
S.conv.N_release = [31 61 121 241];
S.conv.dt_release= [5 2.5 1.25 0.625]*1e-3;
S.conv.pilot_N   = 61;  S.conv.pilot_dt = 2.5e-3;
S.conv.paired_seeds = 1:30;
S.conv.cells = [31 1.25; 61 1.25; 121 1.25; 241 1.25; ...   % mesh at fine step
                121 5; 121 2.5; 121 0.625; ...              % step at fine mesh
                61 2.5; 61 5; 31 2.5; 31 5; ...             % production candidates
                241 0.625];                                 % reference cell
S.conv.cells(:,2) = S.conv.cells(:,2)*1e-3;
S.conv.U_eq_check = [4 6 8 10];                % m/s, mean-equilibrium comparison
% acceptance
S.acc.f1 = 1e-3; S.acc.release_f = 2e-3; S.acc.release_zeta = 0.02;
S.acc.lam = 5e-3; S.acc.outcome = 29; S.acc.site = 29;

% ---- modal evolution with wind (reference case) ----
S.modal_U = 0:0.5:20;

% ---- ensembles at fixed speeds ----
S.ens.sway_factors = [0.8 0.9 1.0 1.1];        % x U50 of the case, stop = false
S.ens.root_grid = [0.7 0.8 0.9 1.0 1.1 1.2 1.35 1.5];  % x U50 of 'root_stem'

% ---- HPC ----
S.max_hours = 70;                              % stop submitting records after this
S.ropts_ref = struct('tail_T', 1.0, 'tail_every', 2, 'hist_every', 4);
end
