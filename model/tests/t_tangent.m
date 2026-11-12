function err = t_tangent(law, support, verbose)
%T_TANGENT Finite-difference check of dR/dq, dR/da and dR/dv.
if nargin < 1, law = 'crossflow'; end
if nargin < 2, support = 'rotational'; end
if nargin < 3, verbose = true; end
P = v35_params(); P.ear.law = law; P.support = support; P.clamped = false;
P.k_theta = 0.644; P.k_M = 0.305; P.k_W = 5.9;
m = v35_build(P, 31);
n = m.n; s = m.s;
r0 = 2*v35_uniform(3, 9, 2*n + 6*n) - 1;
q = zeros(3*n,1);
% a bent, leaning configuration
phi = 0.6*(s/m.L).^1.5; th = pi/2 - phi;
x = cumtrapz(s, cos(th)); y = cumtrapz(s, sin(th));
q(1:3:end) = x + 1e-4*r0(1:n); q(2:3:end) = y + 1e-4*r0(n+1:2*n); q(3:3:end) = th;
q(1:2) = 0;
v = 0.3*r0(2*n+1:5*n); v(1:2) = 0; a = 2*r0(5*n+1:8*n); a(1:2) = 0;
wt = [7.5, -1.2];
f = m.free; nf = m.nf;
[R0, Kq] = v35_residual(m, q, v, a, wt, 0, 0);
[~, KM] = v35_residual(m, q, v, a, wt, 1, 0); Mn = KM - Kq;
[~, KV] = v35_residual(m, q, v, a, wt, 0, 1); Vn = KV - Kq;
hq = 1e-7; ha = 1e-2; hv = 1e-5; Kfd = zeros(nf); Mfd = zeros(nf); Vfd = zeros(nf);
for j = 1:nf
    e = zeros(3*n,1); e(f(j)) = 1;
    d = v35_residual(m, q+hq*e, v, a, wt, 0, 0) - v35_residual(m, q-hq*e, v, a, wt, 0, 0); Kfd(:,j) = d(f)/(2*hq);
    d = v35_residual(m, q, v, a+ha*e, wt, 0, 0) - v35_residual(m, q, v, a-ha*e, wt, 0, 0); Mfd(:,j) = d(f)/(2*ha);
    d = v35_residual(m, q, v+hv*e, a, wt, 0, 0) - v35_residual(m, q, v-hv*e, a, wt, 0, 0); Vfd(:,j) = d(f)/(2*hv);
end
rel = @(A,B) norm(full(A)-B, 'fro')/max(norm(B,'fro'), 1e-300);
err = [rel(Kq, Kfd), rel(Mn, Mfd), rel(Vn, Vfd), norm(full(v35_mass(m,q)) - Mfd,'fro')/norm(Mfd,'fro')];
if verbose
    fprintf('%-10s %-10s dR/dq %.2e  dR/da %.2e  dR/dv %.2e  mass %.2e\n', law, support, err);
end
end
