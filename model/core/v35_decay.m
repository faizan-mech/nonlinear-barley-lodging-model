function D = v35_decay(t, x, xref, ncyc)
%V35_DECAY Frequency and damping ratio from a free-decay signal x(t).
% x is measured about the equilibrium value xref. The frequency is the inverse
% mean interval between up-crossings over the first ncyc cycles (linear
% interpolation); the damping ratio follows from the logarithmic decrement
% of the positive peaks over the same cycles (parabolic peak interpolation).
if nargin < 3 || isempty(xref), xref = 0; end
if nargin < 4 || isempty(ncyc), ncyc = 4; end
t = t(:); x = x(:) - xref;
s = sign(x); iu = find(s(1:end-1) < 0 & s(2:end) >= 0);
tc = t(iu) - x(iu).*(t(iu+1) - t(iu))./(x(iu+1) - x(iu));
tc = tc(1:min(end, ncyc+1));
D = struct('f', NaN, 'zeta', NaN, 'ncycles', max(numel(tc)-1, 0));
if numel(tc) >= 2, D.f = (numel(tc)-1)/(tc(end) - tc(1)); end
ip = find(x(2:end-1) > x(1:end-2) & x(2:end-1) >= x(3:end) & x(2:end-1) > 0) + 1;
ip = ip(t(ip) > tc(1));
ip = ip(1:min(end, ncyc));
if numel(ip) >= 2
    pk = zeros(numel(ip),1);
    for j = 1:numel(ip)
        i = ip(j); y1 = x(i-1); y2 = x(i); y3 = x(i+1);
        den = y1 - 2*y2 + y3;
        if den ~= 0, dd = 0.5*(y1 - y3)/den; else, dd = 0; end
        pk(j) = y2 - 0.25*(y1 - y3)*dd;
    end
    delta = log(pk(1)/pk(end))/(numel(pk)-1);
    D.zeta = delta/sqrt(4*pi^2 + delta^2);
end
end
