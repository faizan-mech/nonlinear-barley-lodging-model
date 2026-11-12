function B = v35_budget(m, q, v, a, wt, k)
%V35_BUDGET Contributions to the bending moment at node k.
%
% The moment at node k is resolved into the moments, about node k, of the
% generalised loads carried by the nodes above it: aerodynamic load, weight
% and inertia of the ear and of the culm, damping, and canopy support.
% B.M(i) are those moments (N m), B.share = B.M/sum(B.M); B.closure is
% sum(B.M) divided by the nodal bending moment used in the failure index.

nd = m.nd; n = m.n; h = m.h;
F = zeros(nd, 8);
names = {'ear_aero','ear_weight','ear_inertia','culm_aero','culm_weight', ...
         'culm_inertia','damping','support'};
tip = m.tip; E = m.ear;
th = q(tip(3)); om = v(tip(3)); ct = cos(th); st = sin(th);
d1 = [-st; ct]; d3 = [ct; st];
% ear
z = zeros(nd,1);
Ra = v35_residual(aero_only_model(m, E), q, v, z, wt, 0, 0);
F(tip,1) = -Ra(tip);
F(tip,2) = [0; -E.swf*E.m*E.g; -E.swc*E.m*E.g*E.c*ct];
mI = E.m*E.si; JI = E.J*E.si; c = E.c; at = a(tip);
F(tip,3) = -[mI*at(1:2) + mI*c*at(3)*d1 - mI*c*om^2*d3; mI*c*(d1.'*at(1:2)) + (mI*c^2 + JI)*at(3)];
% culm drag (as in V35_RESIDUAL)
thn = q(3:3:end); thm = 0.5*(thn(1:end-1) + thn(2:end)); cc = cos(thm); ss = sin(thm);
vx = 0.5*(v(1:3:end-3) + v(4:3:end)); vy = 0.5*(v(2:3:end-3) + v(5:3:end));
relx = m.alpha_e*wt(1) - vx; rely = m.alpha_e*wt(2) - vy;
un = cc.*rely - ss.*relx; fh = 0.5*h.*m.cst.*un.*abs(un);
Fx = zeros(n,1); Fy = zeros(n,1);
Fx(1:end-1) = Fx(1:end-1) - fh.*ss; Fx(2:end) = Fx(2:end) - fh.*ss;
Fy(1:end-1) = Fy(1:end-1) + fh.*cc; Fy(2:end) = Fy(2:end) + fh.*cc;
F(1:3:end,4) = Fx; F(2:3:end,4) = Fy;
F(:,5) = m.Fg;
F(:,6) = -(m.Mc*a);
if ~isempty(m.damp)
    F(m.free,7) = -m.damp.U*(m.damp.cd.*(m.damp.U.'*v(m.free)));
end
W = m.Wn; dt1 = thn(1:end-1) - m.theta_ref; dt2 = thn(2:end) - m.theta_ref;
if m.rot
    r3 = m.k_M*(W(:,1).*dt1 + W(:,2).*dt2); r6 = m.k_M*(W(:,2).*dt1 + W(:,3).*dt2);
    Mt = zeros(n,1); Mt(1:end-1) = Mt(1:end-1) - r3; Mt(2:end) = Mt(2:end) - r6;
    F(3:3:end,8) = Mt;
else
    x1 = q(1:3:end-3); x2 = q(4:3:end);
    r1 = m.k_W*(W(:,1).*x1 + W(:,2).*x2); r4 = m.k_W*(W(:,2).*x1 + W(:,3).*x2);
    Xt = zeros(n,1); Xt(1:end-1) = Xt(1:end-1) - r1; Xt(2:end) = Xt(2:end) - r4;
    F(1:3:end,8) = Xt;
end
X = q(1:3:end); Y = q(2:3:end);
sel = (k+1):n;
lx = X(sel) - X(k); ly = Y(sel) - Y(k);
Mi = zeros(1,8);
for i = 1:8
    Mi(i) = sum(lx.*F(3*sel-1,i) - ly.*F(3*sel-2,i) + F(3*sel,i));
end
Mn = v35_moments(m, q);
B = struct('names', {names}, 'M', Mi, 'total', sum(Mi), 'share', Mi/sum(Mi), ...
    'M_node', Mn(k), 'closure', sum(Mi)/Mn(k), 'node', k, 's', m.s(k));
end

function m2 = aero_only_model(m, E)
% model whose residual is minus the ear aerodynamic load only
m2 = m; m2.Mc = sparse(m.nd, m.nd); m2.Fg = zeros(m.nd,1); m2.cst = 0*m.cst;
m2.k_theta = 0; m2.k_M = 0; m2.k_W = 0; m2.damp = [];
m2.EA = 0*m.EA; m2.GA = 0*m.GA; m2.EI = 0*m.EI;
E.m = 0; E.J = 0; m2.ear = E;
end
