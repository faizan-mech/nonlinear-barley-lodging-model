function [f, Phi, lam] = v35_modes(m, q, wt, k)
%V35_MODES Lowest k natural frequencies of the tangent operator at state q.
% K = dR/dq (including aerodynamic stiffness when wt ~= 0) and M(q) on the
% free dofs. Returns frequencies in Hz (ascending), mass-normalised modes.
if nargin < 3 || isempty(wt), wt = [0 0]; end
if nargin < 4, k = 3; end
z = zeros(m.nd,1);
[~, K] = v35_residual(m, q, z, z, wt, 0, 0);
M = v35_mass(m, q);
[lam, Phi] = v35_eig_low(K, M, k);
f = sqrt(abs(lam))/(2*pi).*sign(lam);
for j = 1:size(Phi,2)
    Phi(:,j) = Phi(:,j)/sqrt(Phi(:,j).'*M*Phi(:,j));
end
end
