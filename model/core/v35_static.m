function [q, ok, it, nr] = v35_static(m, q, wt, opts)
%V35_STATIC Newton solution of the static equilibrium R(q; wt) = 0.
% Full Newton steps, scaled down only when a rotation increment exceeds
% opts.max_dtheta (0.2 rad). The residual norm is not used as a merit
% function: the stiff axial and shear terms make it rise on the first
% step of a geometrically exact rod even when Newton converges.
if nargin < 4, opts = struct(); end
tol = max(getopt(opts, 'tol', 1e-10), m.R_floor);
maxit = getopt(opts, 'max_iter', 30);
dthmax = getopt(opts, 'max_dtheta', 0.2);
f = m.free; z = zeros(m.nd,1);
isth = mod(f, 3) == 0;
ok = false; nr = Inf;
for it = 1:maxit
    [R, K] = v35_residual(m, q, z, z, wt, 0, 0);
    r = R(f); nr = norm(r, inf);
    if nr <= tol, ok = true; return; end
    dq = -(K\r);
    if any(~isfinite(dq)), return; end
    if norm(dq, inf) < 1e-14*m.L && nr < 1e3*tol, ok = true; return; end
    mx = max(abs(dq(isth)));
    if mx > dthmax, dq = dq*(dthmax/mx); end
    q(f) = q(f) + dq;
end
R = v35_residual(m, q, z, z, wt, 0, 0); nr = norm(R(f), inf);
ok = nr <= tol;
end

function v = getopt(o, name, def)
if isfield(o, name), v = o.(name); else, v = def; end
end
