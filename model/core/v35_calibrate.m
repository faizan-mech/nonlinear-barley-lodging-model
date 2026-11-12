function C = v35_calibrate(P, rule, Cref, Ncal)
%V35_CALIBRATE Calibrate the root spring and the canopy support.
%
% rule = 'both' : k_theta so that f1 = P.f_free with no canopy support, then
%                 the canopy modulus so that f1 = P.f_canopy (reference and
%                 alternative constructions of the culm)
% rule = 'kM'   : k_theta = Cref.k_theta*P.k_theta_factor, canopy modulus
%                 recalibrated (window, support-type and k_theta-bound cases)
% rule = 'none' : both from Cref (k_theta times P.k_theta_factor)
%
% Frequencies are those of the tangent operator at the gravity equilibrium,
% computed on the calibration mesh Ncal (default 481 nodes). The canopy
% modulus is k_M (rotational foundation) or k_W (Winkler), per P.support.

if nargin < 2 || isempty(rule), rule = 'both'; end
if nargin < 3, Cref = []; end
if nargin < 4 || isempty(Ncal), Ncal = 481; end
C = struct('rule', rule, 'Ncal', Ncal, 'support', P.support);

if strcmp(rule, 'none')
    C.k_theta = Cref.k_theta*P.k_theta_factor;
    C.k_M = Cref.k_M; C.k_W = Cref.k_W;
    return;
end

% Linear dependence of the tangent on the supports at the upright state:
% K = K0 + k_theta*e3*e3' + k_mod*Kf.
P0 = P; P0.k_theta = 0; P0.k_M = 0; P0.k_W = 0; P0.clamped = false;
m = v35_build(P0, Ncal);
q = zeros(m.nd,1); q(2:3:end) = m.s; q(3:3:end) = pi/2;
m.k_theta = 1e3;                     % any positive value for the gravity solve
[qg, ok] = v35_static(m, q, [0 0]);
if ~ok, error('V35:calibrate','Gravity equilibrium failed.'); end
z = zeros(m.nd,1);
m.k_theta = 0; [~, K0] = v35_residual(m, qg, z, z, [0 0], 0, 0);
E3 = sparse(m.fmap(3), m.fmap(3), 1, m.nf, m.nf);
if strcmpi(P.support, 'rotational'), m.k_M = 1; else, m.k_W = 1; end
[~, K1] = v35_residual(m, qg, z, z, [0 0], 0, 0);
Kf = K1 - K0;
M = v35_mass(m, qg);
f1 = @(kt, km) sqrt(max(v35_eig_low(K0 + kt*E3 + km*Kf, M, 1), 0))/(2*pi);
fk = @(kt, km, k) sqrt(max(v35_eig_low(K0 + kt*E3 + km*Kf, M, k), 0))/(2*pi);

% clamped reference
Pc = P0; Pc.clamped = true; mc = v35_build(Pc, Ncal); mc = v35_prepare(mc);
C.f_clamp = mc.f_grav(1:3).';

if strcmp(rule, 'both')
    C.k_theta = solve_k(@(k) f1(k, 0), P.f_free);
    C.k_theta_interval = [solve_k(@(k) f1(k, 0), P.f_free - P.f_free_band), ...
                          solve_k(@(k) f1(k, 0), P.f_free + P.f_free_band)];
    C.f_free = fk(C.k_theta, 0, 3).';
    C.halving_drop = f1(C.k_theta/2, 0)/P.f_free - 1;
    C.k_theta_ratio_for_20pc_drop = solve_k(@(k) f1(k, 0), 0.8*P.f_free)/C.k_theta;
    C.crook_ennos = P.crook_ennos;
    C.in_crook_ennos = C.k_theta >= P.crook_ennos(1) && C.k_theta <= P.crook_ennos(2);
    C.k_theta_calibrated = C.k_theta;
    C.k_theta = C.k_theta*P.k_theta_factor;
elseif strcmp(rule, 'kM')
    C.k_theta_calibrated = Cref.k_theta;
    C.k_theta = Cref.k_theta*P.k_theta_factor;
    C.f_free = fk(C.k_theta, 0, 3).';
else
    error('V35:calibrate','Unknown rule %s.', rule);
end

if f1(C.k_theta, 0) >= P.f_canopy
    error('V35:calibrate','Free-standing f1 already above the canopy target.');
end
kc = solve_k(@(k) f1(C.k_theta, k), P.f_canopy);
if strcmpi(P.support, 'rotational'), C.k_M = kc; C.k_W = 0; else, C.k_M = 0; C.k_W = kc; end
C.f_canopy = fk(C.k_theta, kc, 3).';
ft = (P.f_canopy-0.10:0.02:P.f_canopy+0.10);
C.canopy_table = [ft; arrayfun(@(x) solve_k(@(k) f1(C.k_theta, k), x), ft)].';
end

function k = solve_k(fun, target)
% Monotone increasing fun(k), k > 0: bracket in log k, then fzero.
g = @(lk) fun(exp(lk)) - target;
a = log(1e-3); b = log(10);
ga = g(a); gb = g(b);
while ga > 0, b = a; gb = ga; a = a - 2; ga = g(a); if a < log(1e-12), error('V35:calibrate','No bracket.'); end, end
while gb < 0, a = b; ga = gb; b = b + 2; gb = g(b); if b > log(1e8), error('V35:calibrate','Target frequency not reachable.'); end, end
lk = fzero(g, [a b], optimset('TolX', 1e-12));
k = exp(lk);
end
