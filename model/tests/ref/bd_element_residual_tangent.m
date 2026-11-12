function [Re,Ke,parts] = bd_element_residual_tangent(qe,h,EA,GA,EI)
%BD_ELEMENT_RESIDUAL_TANGENT Exact 2-node planar Cosserat element.
%
% Unknown order: [x1 y1 theta1 x2 y2 theta2]. One-point axial/shear
% integration avoids locking; curvature is constant and integrated exactly.

r1 = qe(1:2); th1 = qe(3);
r2 = qe(4:5); th2 = qe(6);
rp = (r2-r1)/h;
thm = 0.5*(th1+th2);
d3 = [cos(thm);sin(thm)];
d1 = [-sin(thm);cos(thm)];

epsa = d3.'*rp - 1;
gam  = d1.'*rp;
kappa = (th2-th1)/h;

G = [-1/h 0 0 1/h 0 0; 0 -1/h 0 0 1/h 0];
a = [0;0;0.5;0;0;0.5];
b = [0;0;-1/h;0;0;1/h];

ge = G.'*d3 + gam*a;
gg = G.'*d1 - (epsa+1)*a;
He = (G.'*d1)*a.' + a*(G.'*d1).' - (epsa+1)*(a*a.');
Hg = -(G.'*d3)*a.' - a*(G.'*d3).' - gam*(a*a.');

Re = h*(EA*epsa*ge + GA*gam*gg + EI*kappa*b);
Ke = h*(EA*(ge*ge.' + epsa*He) + ...
        GA*(gg*gg.' + gam*Hg) + EI*(b*b.'));
Ke = 0.5*(Ke+Ke.');

if nargout>2
    parts = struct();
    parts.eps_axial = epsa;
    parts.gamma = gam;
    parts.kappa = kappa;
    parts.energy_axial = 0.5*h*EA*epsa^2;
    parts.energy_shear = 0.5*h*GA*gam^2;
    parts.energy_bending = 0.5*h*EI*kappa^2;
    parts.K_bending = h*EI*(b*b.');
end
end
