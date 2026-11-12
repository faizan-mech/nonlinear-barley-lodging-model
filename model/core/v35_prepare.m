function m = v35_prepare(m)
%V35_PREPARE Gravity equilibrium, natural frequencies and modal damping.
%
% The damping matrix is C = M0*Phi*diag(2*zeta*omega)*Phi'*M0 over the
% first P.n_damped modes of the tangent operator at the gravity equilibrium
% (still air), so each of these modes has damping ratio zeta exactly. It is
% held constant during a record and applied in low-rank form.
P = m.P;
if ~P.clamped && ~(isfinite(m.k_theta) && m.k_theta > 0)
    error('V35:prepare','k_theta is not set (calibrate first).');
end
q = zeros(m.nd,1); q(2:3:end) = m.s; q(3:3:end) = pi/2;
m.damp = [];
[qg, ok] = v35_static(m, q, [0 0]);
if ~ok, error('V35:prepare','Gravity equilibrium failed.'); end
m.q_grav = qg;
nm = min(P.n_damped, m.nf);
[f, Phi] = v35_modes(m, qg, [0 0], max(nm, 5));
if any(f(1:nm) <= 0), error('V35:prepare','Unstable gravity equilibrium (f1 = %.4g Hz).', f(1)); end
M0 = v35_mass(m, qg);
m.f_grav = f;
if P.zeta > 0
    m.damp = struct('U', M0*Phi(:,1:nm), 'cd', 2*P.zeta*2*pi*f(1:nm), ...
        'f', f(1:nm), 'zeta', P.zeta);
end
end
