function [Uq, info] = v35_qs_threshold(m, Iu, Umax)
%V35_QS_THRESHOLD Speed at which the mean-load equilibrium reaches lambda = 1.
% The mean load is that of Ueq(U, Iu*U); Iu = 0 gives the steady-wind
% threshold. Returns Inf if lambda < 1 up to Umax (default 100 m/s).
if nargin < 2, Iu = 0; end
if nargin < 3, Umax = 100; end
lamU = @(U) lam_at(m, v35_ueq(U, Iu*U));
lo = 0; hi = 4;
while lamU(hi) < 1
    lo = hi; hi = 2*hi;
    if hi > Umax, Uq = Inf; info = struct('lo', lo); return; end
end
for k = 1:40
    mid = 0.5*(lo + hi);
    if lamU(mid) >= 1, hi = mid; else, lo = mid; end
    if hi - lo < 1e-4*hi, break; end
end
Uq = 0.5*(lo + hi);
info = struct('lo', lo, 'hi', hi);
end

function L = lam_at(m, Ueq)
[q, inf0] = v35_equilibrium(m, Ueq);
if ~inf0.ok, L = Inf; return; end
[Mn, Mb] = v35_moments(m, q);
L = -inf;
if m.fail.stem, L = max(abs(Mn)./m.BS); end
if m.fail.root, L = max(L, abs(Mb)/m.fail.root_cap); end
end
