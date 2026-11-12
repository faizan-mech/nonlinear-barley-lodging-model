function wt = v35_wind_at(W, t)
%V35_WIND_AT Air velocity [u w] at ear height at time t (ramp included).
% u = Ueq + r(t)*(U + u'(t) - Ueq), w = r(t)*w'(t), with the ramp
% r = (1 - cos(pi t/T_ramp))/2 from the mean-load equilibrium.
if W.U == 0, wt = [0 0]; return; end
x = t/W.dt; i = floor(x); fr = x - i;
i1 = min(i+1, W.N); i2 = min(i+2, W.N);
u = W.u(i1) + fr*(W.u(i2) - W.u(i1));
w = W.w(i1) + fr*(W.w(i2) - W.w(i1));
if t < W.T_ramp, r = 0.5*(1 - cos(pi*t/W.T_ramp)); else, r = 1; end
wt = [W.Ueq + r*(u - W.Ueq), r*w];
end
