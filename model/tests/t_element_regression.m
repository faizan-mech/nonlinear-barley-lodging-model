function err = t_element_regression(verbose)
%T_ELEMENT_REGRESSION Compare the V35 elastic assembly with the
% reference element (bd_element_residual_tangent, in tests/ref).
if nargin < 1, verbose = true; end
P = v35_params(); P.k_theta = 0.6; P.k_M = 0;
m = v35_build(P, 31);
m.cst(:) = 0; m.k_theta = 0; m.Fg(:) = 0; m.ear.m = 0; m.ear.J = 0;
m.ear.kn = 0; m.ear.ka = 0; m.damp = [];
st = rand('state'); rand('state', 7);
n = m.n; s = m.s;
th = pi/2 - 0.8*(s/m.L).^2 + 0.02*(rand(n,1) - 0.5);
q = zeros(3*n,1);
q(1:3:end) = cumtrapz(s, cos(th)).*(1 + 1e-3*(rand(n,1)-0.5));
q(2:3:end) = cumtrapz(s, sin(th)).*(1 + 1e-3*(rand(n,1)-0.5));
q(3:3:end) = th;
rand('state', st);
z = zeros(3*n,1);
[R, K] = v35_residual(m, q, z, z, [0 0], 0, 0);
Rr = zeros(3*n,1); Kr = zeros(3*n);
for e = 1:m.ne
    d = 3*e-2:3*e+3;
    [Re, Ke] = bd_element_residual_tangent(q(d), m.h(e), m.EA(e), m.GA(e), m.EI(e));
    Rr(d) = Rr(d) + Re; Kr(d,d) = Kr(d,d) + Ke;
end
f = m.free;
err = [norm(R(f) - Rr(f))/norm(Rr(f)), norm(full(K) - Kr(f,f), 'fro')/norm(Kr(f,f), 'fro')];
if verbose
    fprintf('element regression vs reference element: residual %.2e, tangent %.2e\n', err);
end
end
