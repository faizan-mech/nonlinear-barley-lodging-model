function s = v35_mesh(P, nNodes)
%V35_MESH Nested non-uniform mesh with nodes on every section edge.
%
% The 31-node mesh has [7 4 4 4 4 7] elements per segment.
% Successive refinements halve element lengths and retain all coarser nodes.

L  = P.L;
bp = [0, 0.20, 0.315, P.s_taper*L, 0.515, 0.63, L];
if any(diff(bp) <= 0), error('V35:mesh','Breakpoints are not increasing.'); end
base = [7 4 4 4 4 7];
lev = log2((nNodes-1)/sum(base));
if abs(lev-round(lev)) > 1e-12 || lev < 0
    error('V35:mesh','nNodes must be 30*2^k+1 (31, 61, 121, 241, 481); got %d.', nNodes);
end
ne = base*2^round(lev);
s = zeros(sum(ne)+1,1); k = 1;
for j = 1:numel(ne)
    x = linspace(bp(j), bp(j+1), ne(j)+1).';
    s(k:k+ne(j)) = x; k = k+ne(j);
end
s(end) = L;
end
