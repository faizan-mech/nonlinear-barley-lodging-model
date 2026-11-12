function S = v35_sections(P, s, E, rho_s)
%V35_SECTIONS Section properties at arc lengths s (column vector, m).
%
% S = v35_sections(P,s)            geometry and strength only
% S = v35_sections(P,s,E,rho_s)    also stiffness and mass per unit length
%
% Geometry: outer radius uniform to s_taper*L, then linear to taper*R_mid
% at the ear (radius only; wall constant above the middle section). Below
% 0.315 m the wall thickens in proportion to Berry's strength layout.

s = s(:); L = P.L;
Rref = P.R_mid*(1 + (P.taper-1)*max(0,(s-P.s_taper*L)/((1-P.s_taper)*L)));
tref = P.t_mid*wall_factor(P, s);
tref = min(tref, 0.98*Rref);

switch lower(P.wall_keep)
    case 'outer'
        R = P.scale_Do*Rref;
        t = min(P.scale_t*tref, 0.98*R);
    case 'inner'
        if P.scale_Do ~= 1
            error('V35:sections','wall_keep=inner is defined for scale_Do = 1 only.');
        end
        Ri = Rref - tref;
        t  = P.scale_t*tref;
        R  = Ri + t;
    otherwise
        error('V35:sections','Unknown wall_keep %s.', P.wall_keep);
end
Ri = R - t;
S.R = R; S.t = t; S.D = 2*R;
S.A = pi*(R.^2 - Ri.^2);
S.I = pi/4*(R.^4 - Ri.^4);
S.Z = S.I./R;

% Reference section (for BS proportional to Z)
Riref = Rref - tref;
Iref  = pi/4*(Rref.^4 - Riref.^4);
Zref  = Iref./Rref;

BS = P.scale_BS*interp1(P.BS_s(:), P.BS_v(:), min(max(s,P.BS_s(1)),P.BS_s(end)), 'linear');
if P.BS_scales_with_Z
    BS = BS.*S.Z./Zref;
end
S.BS = BS;

w = 0.5*(tanh((s - P.win_lo*L)/(P.win_d*L)) - tanh((s - P.win_hi*L)/(P.win_d*L)));
w(s >= P.win_cut*L) = 0;
S.window = w;
S.alpha  = exp(-P.a_canopy*(1 - s/L));

if nargin >= 3 && ~isempty(E)
    I_bend = S.I;
    if P.EI_follows_strength
        % Upper bound: EI follows the strength of the uniform-wall section.
        Ru = R; Riu = Ru - min(P.t_mid, 0.98*Ru);
        Iu = pi/4*(Ru.^4 - Riu.^4);
        I_bend = Iu.*wall_factor(P, s);
    end
    G = E/(2*(1+P.nu));
    S.EI = E*I_bend;
    S.EA = E*S.A;
    S.GA = P.k_shear*G*S.A;
end
if nargin >= 4 && ~isempty(rho_s)
    sg = P.sig_leaf*L;
    b1 = exp(-0.5*((s - P.s_flag*L)/sg).^2);
    b2 = exp(-0.5*((s - P.s_next*L)/sg).^2);
    S.mu = rho_s*S.A + P.m_shoot*P.frac_flag*b1/bump_norm(P.s_flag*L, sg, L) ...
                     + P.m_shoot*P.frac_next*b2/bump_norm(P.s_next*L, sg, L);
    S.mu_stalk = rho_s*S.A;
    S.jrot = rho_s*S.I;
end
end

function f = wall_factor(P, s)
% Wall-thickness factor t/t_mid of the reference construction.
f = ones(size(s));
if strcmpi(P.lower_wall,'strength')
    lo = s < P.wall_s(end);
    bs = interp1(P.wall_s(:), P.wall_BS(:), min(max(s(lo),P.wall_s(1)),P.wall_s(end)), 'linear');
    f(lo) = bs/P.BS_mid;
end
end

function z = bump_norm(s0, sg, L)
% Integral over [0,L] of exp(-0.5((s-s0)/sg)^2).
z = sg*sqrt(pi/2)*(erf((L-s0)/(sqrt(2)*sg)) - erf((0-s0)/(sqrt(2)*sg)));
end
