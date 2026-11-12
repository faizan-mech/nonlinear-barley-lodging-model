function [R, K] = v35_residual(m, q, v, a, wt, cM, cV)
%V35_RESIDUAL Dynamic residual and effective tangent of the shoot.
%
%   R = M(q) a + h(q,v) + C v + Rint(q) + Rsupport(q) - Fgrav(q) - Faero(q,v,t)
%   h(q,v) contains the velocity-quadratic inertia of the offset rigid head.
%   K = dR/dq + cM*dR/da + cV*dR/dv      (free dofs, sparse; the low-rank
%                                         modal damping is excluded, see V35_SOLVE)
%
% q, v, a : full dof vectors (3n); wt = [u w] air velocity at ear height
% (ramp included). The culm sees alpha(s)*wt at its undeformed height.
% For statics use v = a = 0, cM = cV = 0.
%
% One-point element quadrature; assembly uses m.Sr and m.Sk from V35_BUILD.

h = m.h; n = m.n; ne = m.ne;
th = q(3:3:end);
dx = q(4:3:end) - q(1:3:end-3); dy = q(5:3:end) - q(2:3:end-3);
dth = th(2:end) - th(1:end-1);
thm = 0.5*(th(1:end-1) + th(2:end));
c = cos(thm); s = sin(thm);
ep  = (c.*dx + s.*dy)./h - 1;          % axial strain
gam = (c.*dy - s.*dx)./h;              % shear strain
EA = m.EA; GA = m.GA;
Nf = EA.*ep; Qf = GA.*gam; Mb = m.EI.*dth./h;

% culm crossflow drag, one-point at the element centre
vx = 0.5*(v(1:3:end-3) + v(4:3:end)); vy = 0.5*(v(2:3:end-3) + v(5:3:end));
relx = m.alpha_e*wt(1) - vx; rely = m.alpha_e*wt(2) - vy;
un = c.*rely - s.*relx; aun = abs(un);
fh = 0.5*h.*m.cst.*un.*aun;           % (h/2)*f, f along d1 = (-s, c)

% supports
dt1 = th(1:end-1) - m.theta_ref; dt2 = th(2:end) - m.theta_ref;
W = m.Wn;
if m.rot
    kM = m.k_M;
    r3 = kM*(W(:,1).*dt1 + W(:,2).*dt2); r6 = kM*(W(:,2).*dt1 + W(:,3).*dt2);
    r1 = 0; r4 = 0;
else
    kW = m.k_W; x1 = q(1:3:end-3); x2 = q(4:3:end);
    r1 = kW*(W(:,1).*x1 + W(:,2).*x2); r4 = kW*(W(:,2).*x1 + W(:,3).*x2);
    r3 = 0; r6 = 0;
end
% element residual: [ex ey e36 eb rx ry e1 e3 e4 e6] patterns
Re = [Nf.*c - Qf.*s, Nf.*s + Qf.*c, 0.5*h.*(Nf.*gam - Qf.*(1+ep)), Mb, ...
      fh.*s, -fh.*c, r1 + 0*h, r3 + 0*h, r4 + 0*h, r6 + 0*h]*m.Sr;
Rn = zeros(n,3);
Rn(1:ne,:) = Re(:,1:3);
Rn(2:n,:) = Rn(2:n,:) + Re(:,4:6);
R = reshape(Rn.', [], 1);
R = R + m.Mc*a - m.Fg;
if ~m.clamped
    R(3) = R(3) + m.k_theta*(th(1) - m.theta_ref);
end

% ear (exact rigid body)
[Rt, Kt] = ear_terms(m.ear, q(m.tip), v(m.tip), a(m.tip), wt, cM, cV, nargout > 1);
R(m.tip) = R(m.tip) + Rt;

% modal damping (low rank)
if ~isempty(m.damp)
    vf = v(m.free);
    R(m.free) = R(m.free) + m.damp.U*(m.damp.cd.*(m.damp.U.'*vf));
end
if nargout < 2, return; end

% element tangent coefficients (see V35_BUILD for the 18 patterns)
EAh = EA./h; GAh = GA./h; cc = c.*c; ss = s.*s; cs = c.*s;
dEG = EA - GA; t1 = GA.*(1+ep) - Nf;
ut = c.*relx + s.*rely;
gp = 2*m.cst.*aun;                     % d|f|/d un
gx = gp.*ut.*s - fh./(0.5*h).*c;       % dF/dthm, F = f d1
gy = -gp.*ut.*c - fh./(0.5*h).*s;
q4 = 0.25*h;
kd = cV*0.5*h.*m.cst.*aun;
z = zeros(ne,1);
if m.rot
    sup = [m.k_M*W, z, z, z];
else
    sup = [z, z, z, m.k_W*W];
end
Phi = [EAh.*cc + GAh.*ss, EAh.*ss + GAh.*cc, (EAh - GAh).*cs, ...
       0.5*(c.*dEG.*gam + s.*t1), 0.5*(s.*dEG.*gam - c.*t1), ...
       0.25*h.*(dEG.*gam.^2 + (1+ep).*t1), m.EI./h, ...
       -q4.*gx, -q4.*gy, kd.*ss, kd.*cc, -kd.*cs, sup];
Ke = Phi*m.Sk;
if cM ~= 0, Ke = Ke + cM*m.McV; end
if m.clamped, kb = 0; else, kb = m.k_theta; end
V = [Ke(:); Kt(:); kb];
K = sparse(m.K_I, m.K_J, V(m.K_keep), m.nf, m.nf);
end

function [Rt, Kt] = ear_terms(E, qt, vt, at, wt, cM, cV, wantK)
% Residual contribution and effective tangent (3x3) of the rigid ear.
th = qt(3); om = vt(3);
ct = cos(th); st = sin(th);
d3 = [ct; st]; d1 = [-st; ct];
mI = E.m*E.si; JI = E.J*E.si; c = E.c;
axy = at(1:2); ath = at(3);
% inertia: M(th) a + centripetal
Rt = [mI*axy + mI*c*ath*d1 - mI*c*om^2*d3; mI*c*(d1.'*axy) + (mI*c^2 + JI)*ath];
% weight: force at the tip and couple about the neck
mg = E.m*E.g;
Rt(2) = Rt(2) + E.swf*mg;
Rt(3) = Rt(3) + E.swc*mg*c*ct;
% aerodynamic load
if E.sa ~= 0
    [Fa, dFq, dFv] = ear_aero(E, th, vt, wt);
    Rt = Rt - E.sa*Fa;
else
    dFq = zeros(3); dFv = zeros(3);
end
Kt = zeros(3);
if ~wantK, return; end
Kq = zeros(3);
Kq(1:2,3) = -mI*c*ath*d3 - mI*c*om^2*d1;
Kq(3,3)   = -mI*c*(d3.'*axy) - E.swc*mg*c*st;
Kq = Kq - E.sa*dFq;
Kv = zeros(3);
Kv(1:2,3) = -2*mI*c*om*d3;
Kv = Kv - E.sa*dFv;
Mt = [mI 0 -mI*c*st; 0 mI mI*c*ct; -mI*c*st mI*c*ct mI*c^2 + JI];
Kt = Kq + cM*Mt + cV*Kv;
end

function [F, dFq, dFv] = ear_aero(E, th, vt, wt)
% Generalised aerodynamic force [Fx; Fy; M] on the tip dofs and its
% derivatives with respect to (x,y,theta) and (vx,vy,omega).
ct = cos(th); st = sin(th); om = vt(3);
d3 = [ct; st]; d1 = [-st; ct];
dFq = zeros(3); dFv = zeros(3);
switch E.law
    case 'crossflow'
        % normal force at the centre of pressure, axial force through the neck
        rel = [wt(1); wt(2)] - vt(1:2) - om*E.rcp*d1;
        un = d1.'*rel; ua = d3.'*rel;
        fn = E.kn*un*abs(un); fa = E.ka*ua*abs(ua);
        F = [fn*d1 + fa*d3; E.rcp*fn];
        gn = 2*E.kn*abs(un); ga = 2*E.ka*abs(ua);
        dun_dth = -ua; dua_dth = un + om*E.rcp;
        dFq(1:2,3) = gn*dun_dth*d1 - fn*d3 + ga*dua_dth*d3 + fa*d1;
        dFq(3,3) = E.rcp*gn*dun_dth;
        % dun/dv = -d1', dun/dom = -rcp; dua/dv = -d3', dua/dom = 0
        dFv(1:2,1:2) = -gn*(d1*d1.') - ga*(d3*d3.');
        dFv(1:2,3) = -gn*E.rcp*d1;
        dFv(3,1:2) = -E.rcp*gn*d1.';
        dFv(3,3) = -E.rcp^2*gn;
    case {'projection','area','upright'}
        % horizontal drag at the centre of pressure; area and arm follow the
        % projection (bound), only the area follows it with the arm held at
        % its upright value ('area'), or both stay upright ('upright')
        ur = wt(1) - vt(1) + om*E.rcp*st;         % horizontal relative velocity at CP
        if strcmp(E.law,'upright')
            A = E.kn; dA = 0;
        else
            A  = E.kn*abs(st) + E.ka*abs(ct);
            dA = E.kn*sign(st)*ct - E.ka*sign(ct)*st;
        end
        if strcmp(E.law,'projection')
            arm = E.rcp*st; darm = E.rcp*ct;
        else
            arm = E.rcp; darm = 0;
        end
        Fx = A*ur*abs(ur); g = 2*A*abs(ur);
        F = [Fx; 0; -arm*Fx];
        dFx_dth = dA*ur*abs(ur) + g*om*E.rcp*ct;
        dFq(1,3) = dFx_dth;
        dFq(3,3) = -darm*Fx - arm*dFx_dth;
        dFv(1,1) = -g; dFv(1,3) = g*E.rcp*st;
        dFv(3,1) = arm*g; dFv(3,3) = -arm*g*E.rcp*st;
end
end
