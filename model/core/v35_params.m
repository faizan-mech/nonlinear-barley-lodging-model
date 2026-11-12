function P = v35_params()
%V35_PARAMS Reference physical inputs of the V35 barley shoot model.
%
% Sensitivity cases apply overrides through V35_APPLY_VARIANT.
% Lengths in m, masses in kg, forces in N.
%
% Sources:
%   Berry et al. (2006)            culm length, middle-20-cm EI, section
%                                  strengths, free-standing frequency 0.97 Hz
%   Neenan & Spencer-Smith (1975)  outer-radius taper ratio 0.57
%   Ren (2014)                     internode proportions 36/22/42 %
%   Finnigan (1979); AHDB          mass fractions of stalk and leaves; ear mass
%   Susko (2019)                   supported frequency 1.36 Hz
%   Finnigan & Mulhearn (1978)     exponential canopy profile, a = 2.6
%   Baker et al. (1998)            root anchorage M_R = k3*tau_s*d^3

P = struct();
P.version = 'V35';

%% Culm geometry
P.L          = 0.83;        % culm length, cv. Pearl (Berry et al. 2006)
P.R_mid      = 1.8857e-3;   % outer radius, uniform up to the taper start
P.t_mid      = 0.39e-3;     % wall thickness of the middle section
P.s_taper    = 0.5;         % taper start (fraction of L)
P.taper      = 0.57;        % R(L)/R_mid, radius only; wall constant (Neenan Table 1)
P.s_ped      = 0.58;        % peduncle base (fraction of L), Ren (2014) 36/22/42 %
P.nu         = 0.30;        % Poisson ratio
P.k_shear    = 0.55;        % representative hollow-section shear factor (Cowper 1966)

% Lower wall: t = t_mid*BS/BS_mid for s < 0.315 m.
P.lower_wall = 'strength';  % 'strength' (reference) | 'uniform'
P.wall_s     = [0 0.20 0.315];        % m, Berry layout below the middle section
P.wall_BS    = [0.075 0.075 0.042];   % N m; construction profile, fixed in capacity variants

% Lower-culm EI upper bound (sensitivity): EI(s) = EI_uniform(s)*BS/BS_mid
% for s < 0.315 m while geometry and mass stay those of the reference.
P.EI_follows_strength = false;

% Reference E matches the harmonic-mean EI over the middle 20 cm.
% Geometry variants retain this reference value.
P.EI_mid20   = 0.0246;
P.mid20      = [0.315 0.515];   % m
P.E          = [];          % [] -> derived by V35_BUILD

% Geometry perturbations (sensitivity cases); 1 = reference
P.scale_Do         = 1;     % outer radius factor
P.scale_t          = 1;     % wall thickness factor
P.wall_keep        = 'outer';  % with scale_t: 'outer' radius fixed or 'inner' radius fixed
P.scale_E          = 1;

%% Section strength BS(s), Berry et al. (2006) layout for the 0.83 m stem
P.BS_s   = [0 0.20 0.315 0.515 0.63 0.83];   % m
P.BS_v   = [0.075 0.075 0.042 0.042 0.021 0.021];  % N m
P.BS_mid = 0.042;
P.scale_BS = 1;             % strength factor (sensitivity)
P.BS_scales_with_Z = false; % true: BS multiplied by Z/Z_ref (diameter cases)

%% Culm and leaf mass
P.m_shoot    = 3.11e-3;     % culm + sheaths + leaves
P.frac_stalk = 48/62.5;     % stalk, follows its own cross-section area
P.frac_flag  = 8/62.5;      % flag leaf, Gaussian at 0.64 L
P.frac_next  = 6.5/62.5;    % next leaf, Gaussian at 0.45 L
P.s_flag     = 0.64;
P.s_next     = 0.45;
P.sig_leaf   = 0.04;        % Gaussian width (fraction of L)

%% Ear: exact rigid body attached at the culm tip
P.ear.m      = 1.875e-3;    % kg
P.ear.J_G    = 1.1342711191582338e-6;  % kg m^2 about the centre of gravity
P.ear.c      = 41.2e-3;     % centre of gravity along the ear axis beyond the neck
P.ear.r_cp   = 63.2e-3;     % centre of pressure along the ear axis
P.ear.A_up   = 1.605e-3;    % projected area normal to the axis (upright, broadside)
P.ear.A_hor  = 9.39e-5;     % projected area along the axis (horizontal ear, end-on)
P.ear.C_D    = 0.866;
P.ear.law    = 'crossflow'; % 'crossflow' (reference) | 'projection' | 'area' | 'upright'
% Mass pathways (1 = physical). Decomposition cases scale them independently.
P.ear.s_inertia      = 1;   % translational + rotary inertia (incl. coupling)
P.ear.s_weight_force = 1;   % weight carried at the culm tip
P.ear.s_weight_couple= 1;   % couple of the weight about the neck, -m g c cos(theta)
P.ear.s_aero         = 1;

%% Air and culm drag
P.rho_air  = 1.225;
P.C_D_culm = 1.0;           % crossflow drag on the culm (outer diameter)
P.g        = 9.81;
P.a_canopy = 2.6;           % alpha(s) = exp[-a(1-s/L)] at undeformed height

%% Supports
P.f_free   = 0.97;          % free-standing first frequency (k_M = 0), Hz
P.f_free_band = 0.025;      % +/- band for the k_theta interval
P.f_canopy = 1.36;          % in-canopy first frequency, Hz
P.k_theta  = [];            % N m/rad, [] -> calibrated
P.k_M      = [];            % N/rad (rotational foundation modulus), [] -> calibrated
P.k_W      = 0;             % N/m^2 (lateral Winkler modulus), variant only
P.support  = 'rotational';  % 'rotational' (reference) | 'winkler'
P.win_lo   = 0.30;          % lower transition centre of the support window (fraction of L)
P.win_hi   = 0.625;         % upper transition centre of the support window (fraction of L)
P.win_d    = 0.03;          % tanh transition width (fraction of L)
P.win_cut  = 0.75;          % hard cut-off (fraction of L)
P.clamped  = false;         % rigid base (calibration check only)
P.k_theta_factor = 1;       % applied after calibration (k_theta bound cases)
P.crook_ennos = [0.57 0.78];% N m/rad, whole-plant wheat secant stiffnesses at 30 degrees

%% Damping
P.zeta     = 0.08;          % modal damping ratio
P.n_damped = 10;            % number of damped modes (gravity equilibrium)

%% Failure
P.fail.stem      = true;
P.fail.root      = true;
P.fail.k3        = 0.58;    % Baker et al. (1998)
P.fail.d_root    = 44e-3;   % root-plate spread, m
P.fail.tau_s     = 4.0e3;   % soil shear strength (field), Pa
P.fail.N_shoots  = 4.4;     % shoots per plant
P.fail.root_share= false;   % true: capacity M_R/N (root runs); false: M_R
P.fail.root_cap  = [];      % N m; [] -> k3*tau_s*d^3 (/N if root_share)

%% Inflow
P.wind.spectrum = 'vonkarman';  % 'vonkarman' (reference) | 'kaimal'
P.wind.Iu       = 0.40;
P.wind.LuL      = 1.0;          % Lu/L
P.wind.vertical = false;        % correlated w' (sensitivity)
P.wind.sw_su    = 0.55;         % sigma_w/sigma_u
P.wind.r_uw     = -0.45;
P.wind.LwL      = 0.35;         % Lw/L (0.25-0.5 h)
P.wind.T        = 60;           % s, record length
P.wind.T_ramp   = 2;            % s, ramp from the mean-load equilibrium
P.wind.dt       = 0.625e-3;     % s, synthesis step (independent of the solver step)
end
