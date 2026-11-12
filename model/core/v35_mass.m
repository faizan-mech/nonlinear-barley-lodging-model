function M = v35_mass(m, q)
%V35_MASS Mass matrix on the free dofs at configuration q (culm + exact ear).
E = m.ear; th = q(3*m.n);
mI = E.m*E.si; JI = E.J*E.si; c = E.c;
Mt = [mI 0 -mI*c*sin(th); 0 mI mI*c*cos(th); -mI*c*sin(th) mI*c*cos(th) mI*c^2 + JI];
Me = sparse(m.nd, m.nd);
Me(m.tip, m.tip) = Mt;
M = m.Mc_ff + Me(m.free, m.free);
end
