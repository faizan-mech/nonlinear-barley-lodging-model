function [Mn, Mb] = v35_moments(m, q)
%V35_MOMENTS Bending moment at the nodes and the anchorage moment.
% Element moments EI*kappa (constant per element) are averaged to interior
% nodes; the base node carries the root-spring moment, the tip node the
% moment of the last element.
th = q(3:3:end);
Me = m.EI.*diff(th)./m.h;
if m.P.clamped
    Mb = Me(1);
else
    Mb = m.k_theta*(th(1) - m.theta_ref);
end
Mn = [Mb; 0.5*(Me(1:end-1) + Me(2:end)); Me(end)];
end
