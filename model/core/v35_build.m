function m = v35_build(P, nNodes)
%V35_BUILD Finite-element model of the barley shoot for parameters P.
%
% Planar Simo-Reissner rod, two-node elements with one-point (reduced)
% quadrature for axial, shear and bending strain. Consistent mass, gravity
% and foundation matrices are integrated with 5-point Gauss quadrature over
% the exact section fields. The ear is an exact rigid body at the tip.
% Dofs per node: [x y theta], theta measured from horizontal (upright pi/2).

if nargin < 2, nNodes = 61; end
m = struct();
m.P = P;
L = P.L;

%% Modulus and stalk density from the reference construction
Pref = P; Pref.scale_Do = 1; Pref.scale_t = 1; Pref.wall_keep = 'outer';
Pref.EI_follows_strength = false;
if isempty(P.E)
    sm = linspace(P.mid20(1), P.mid20(2), 4001).';
    Sm = v35_sections(Pref, sm);
    E0 = P.EI_mid20*trapz(sm, 1./Sm.I)/(sm(end)-sm(1));
else
    E0 = P.E;
end
E = E0*P.scale_E;
sf = linspace(0, L, 20001).';
Sf = v35_sections(Pref, sf);
rho_s = P.m_shoot*P.frac_stalk/trapz(sf, Sf.A);
m.E = E; m.E_ref = E0; m.rho_s = rho_s;

%% Mesh and quadrature
s = v35_mesh(P, nNodes); n = numel(s); ne = n-1;
h = diff(s);
m.s = s; m.n = n; m.ne = ne; m.h = h; m.L = L; m.nd = 3*n;
[xg, wg] = gauss5();
sg = s(1:end-1) + 0.5*h.*(1 + xg);          % ne x 5
N1 = 0.5*(1 - xg); N2 = 0.5*(1 + xg);       % 1 x 5
Sg = v35_sections(P, sg(:), E, rho_s);
rs = @(f) reshape(f, ne, 5);
avg = @(f) 0.5*(rs(f)*wg.');                % element mean
int = @(f,a,b) 0.5*h.*(rs(f).*(a.*b))*wg.'; % int f Na Nb ds

m.EA = avg(Sg.EA); m.GA = avg(Sg.GA); m.EI = avg(Sg.EI);
Mt = [int(Sg.mu,N1,N1), int(Sg.mu,N1,N2), int(Sg.mu,N2,N2)];
Mr = [int(Sg.jrot,N1,N1), int(Sg.jrot,N1,N2), int(Sg.jrot,N2,N2)];
Wn = [int(Sg.window,N1,N1), int(Sg.window,N1,N2), int(Sg.window,N2,N2)];
fg1 = int(Sg.mu, N1, 1); fg2 = int(Sg.mu, N2, 1);     % int mu N_i ds

smid = s(1:end-1) + 0.5*h;
Smid = v35_sections(P, smid);
m.cst = 0.5*P.rho_air*P.C_D_culm*avg(Sg.D);   % culm crossflow coefficient per length
m.alpha_e = Smid.alpha;                        % canopy profile at element centres
m.smid = smid;

Sn = v35_sections(P, s, E, rho_s);
m.BS = Sn.BS; m.R = Sn.R; m.t = Sn.t; m.D = Sn.D; m.window = Sn.window;
m.EI_node = Sn.EI; m.mu_node = Sn.mu;
m.mass_culm = sum(Mt(:,1) + 2*Mt(:,2) + Mt(:,3));
m.mass_stalk = trapz(sf, v35_sections(P, sf, E, rho_s).mu_stalk);
m.centroid = sum(fg1.*s(1:end-1) + fg2.*s(2:end))/sum(fg1 + fg2);

%% Dof bookkeeping
nd = 3*n;
if P.clamped, fixed = [1 2 3]; else, fixed = [1 2]; end
free = setdiff((1:nd).', fixed(:));
fmap = zeros(nd,1); fmap(free) = 1:numel(free);
m.fixed = fixed(:); m.free = free; m.nf = numel(free); m.fmap = fmap;
m.tip = [3*n-2; 3*n-1; 3*n];
m.theta_ref = pi/2;

e = (1:ne).';
m.edof = [3*e-2, 3*e-1, 3*e, 3*e+1, 3*e+2, 3*e+3];      % ne x 6
Ir = repmat(1:6, 1, 6); Ic = kron(1:6, ones(1,6));
m.Ir = Ir; m.Ic = Ic;
I36 = m.edof(:, Ir); J36 = m.edof(:, Ic);
It = repmat(m.tip, 1, 3); Jt = repmat(m.tip.', 3, 1);
Iall = [I36(:); It(:); 3]; Jall = [J36(:); Jt(:); 3];
keep = fmap(Iall) > 0 & fmap(Jall) > 0;
m.K_keep = keep; m.K_I = fmap(Iall(keep)); m.K_J = fmap(Jall(keep));
m.n36 = 36*ne;

%% Constant element patterns (see V35_RESIDUAL)
ex = [-1 0 0 1 0 0]; ey = [0 -1 0 0 1 0]; e36 = [0 0 1 0 0 1]; eb = [0 0 -1 0 0 1];
rx = [1 0 0 1 0 0]; ry = [0 1 0 0 1 0];
u1 = [1 0 0 0 0 0]; u3 = [0 0 1 0 0 0]; u4 = [0 0 0 1 0 0]; u6 = [0 0 0 0 0 1];
pat = @(r, c) kron(c, r);            % vec of r'*c, column-major
sym = @(r, c) pat(r, c) + pat(c, r);
m.Sr = [ex; ey; e36; eb; rx; ry; u1; u3; u4; u6];
m.Sk = [pat(ex,ex); pat(ey,ey); sym(ex,ey); sym(ex,e36); sym(ey,e36); ...
        pat(e36,e36); pat(eb,eb); pat(rx,e36); pat(ry,e36); ...
        pat(rx,rx); pat(ry,ry); sym(rx,ry); ...
        pat(u3,u3); sym(u3,u6); pat(u6,u6); pat(u1,u1); sym(u1,u4); pat(u4,u4)];
m.clamped = logical(P.clamped);

%% Constant culm mass (element values ne x 36 and global sparse)
McV = zeros(ne, 36);
pos = @(r,c) (c-1)*6 + r;
for d = [1 2]    % x and y translational
    McV(:,pos(d,d))     = Mt(:,1); McV(:,pos(d,d+3)) = Mt(:,2);
    McV(:,pos(d+3,d))   = Mt(:,2); McV(:,pos(d+3,d+3)) = Mt(:,3);
end
McV(:,pos(3,3)) = Mr(:,1); McV(:,pos(3,6)) = Mr(:,2);
McV(:,pos(6,3)) = Mr(:,2); McV(:,pos(6,6)) = Mr(:,3);
m.McV = McV;
m.Mc = sparse(I36(:), J36(:), McV(:), nd, nd);
m.Mc_ff = m.Mc(free, free);

%% Gravity of the culm (constant)
Fg = zeros(nd,1);
Fg(3*(1:n)-1) = -P.g*([fg1; 0] + [0; fg2]);
m.Fg = Fg;

%% Foundation (window integrals; modulus applied at run time)
m.Wn = Wn;                       % [W11 W12 W22] per element
m.support = lower(P.support);
m.k_theta = P.k_theta; m.k_M = P.k_M; m.k_W = P.k_W;
if isempty(m.k_theta), m.k_theta = NaN; end
if isempty(m.k_M), m.k_M = 0; end
if isempty(m.k_W), m.k_W = 0; end
if ~any(strcmp(m.support, {'rotational','winkler'}))
    error('V35:build','Unknown support %s.', P.support);
end
m.rot = strcmp(m.support, 'rotational');

%% Ear
ea = P.ear;
m.ear = struct('m', ea.m, 'J', ea.J_G, 'c', ea.c, 'rcp', ea.r_cp, ...
    'kn', 0.5*P.rho_air*ea.C_D*ea.A_up, 'ka', 0.5*P.rho_air*ea.C_D*ea.A_hor, ...
    'law', lower(ea.law), 'si', ea.s_inertia, 'swf', ea.s_weight_force, ...
    'swc', ea.s_weight_couple, 'sa', ea.s_aero, 'g', P.g);
if ~any(strcmp(m.ear.law, {'crossflow','projection','area','upright'}))
    error('V35:build','Unknown ear law %s.', ea.law);
end

%% Failure
m.fail = P.fail;
if isempty(P.fail.root_cap)
    MR = P.fail.k3*P.fail.tau_s*P.fail.d_root^3;
    if P.fail.root_share, MR = MR/P.fail.N_shoots; end
    m.fail.root_cap = MR;
end

%% Residual round-off floor (axial stiffness times machine precision)
m.R_floor = 10*eps*max(m.EA./h)*L;

%% Damping (set by V35_PREPARE)
m.damp = [];
end

function [x, w] = gauss5()
x = [-0.9061798459386640, -0.5384693101056831, 0, 0.5384693101056831, 0.9061798459386640];
w = [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891];
end
