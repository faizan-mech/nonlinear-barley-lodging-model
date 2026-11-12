function out = v35_newmark(m, q0, W, opts)
%V35_NEWMARK Nonlinear transient by average-acceleration Newmark.
%
% out = v35_newmark(m, q0, W, opts) integrates from the state q0 at rest
% (the mean-wind equilibrium) through the record W (see V35_WIND; W.U = 0
% for still air). Newton iterations use the consistent tangent; the modal
% damping enters through the Woodbury identity (V35_SOLVE).
%
% opts.dt        solver step (s)
% opts.T         end time (default W.T)
% opts.stop      stop at first passage of lambda = 1 (default true)
% opts.tol       residual tolerance (default max(1e-9, round-off floor))
% opts.max_iter  Newton iterations per step (12)
% opts.max_cut   recursive step halvings when Newton fails (5)
% opts.ev_tol    accept the event when 1 <= lambda <= 1+ev_tol (2e-3)
% opts.ev_depth  maximum sub-step halvings in event localisation (14)
% opts.ev_guard  accept a localised event with lambda <= guard at max depth (1.05)
% opts.stat_t0   start of the statistics window (default W.T_ramp)
% opts.hist_every  store a decimated history every k steps (0 = none)
% opts.tail_T, opts.tail_every  keep the last tail_T seconds of states
%
% First passage: lambda = max( max_s |M|/BS , |M_root|/M_R ) over the
% active criteria (m.fail.stem, m.fail.root).

o = defaults(opts, W, m);
dt = o.dt; nsteps = round(o.T/dt);
f = m.free; n = m.n;
BSinv = 1./m.BS; capinv = 1/m.fail.root_cap;
useS = m.fail.stem; useR = m.fail.root;

% initial acceleration
q = q0; v = zeros(m.nd,1); a = zeros(m.nd,1);
wt0 = v35_wind_at(W, 0);
R0 = v35_residual(m, q, v, a, wt0, 0, 0);
a(f) = -(v35_mass(m, q)\R0(f));

% statistics
S = struct('cnt', 0, 'sx', zeros(n,1), 'sxx', zeros(n,1), 'sy', zeros(n,1), ...
    'xmin', inf(n,1), 'xmax', -inf(n,1), 'sl', zeros(n,1), 'sll', zeros(n,1), ...
    'sth', 0, 'sthth', 0, 'smb', 0, 'smbmb', 0);
[lamS, lamR, Mb] = lambda_of(m, q, BSinv, capinv);
env = lamS;
[lsmax, ism] = max(lamS); tsm = 0; lrmax = lamR; trm = 0; mbmax = abs(Mb);
Lprev = crit(lamS, lamR, useS, useR);
tcross = NaN;
nh = 0;
if o.hist_every > 0
    H = zeros(floor(nsteps/o.hist_every)+2, 8, 'single');
    nh = 1; H(1,:) = single([0, q(3*n-2), q(3*n-1), q(3*n), lsmax, ism, Mb, wt0(1)]);
else
    H = zeros(0, 8, 'single');
end
tailOn = o.tail_T > 0;
if tailOn
    nb = ceil(o.tail_T/(dt*o.tail_every)) + 1;
    TB = struct('t', nan(1,nb), 'q', zeros(m.nd,nb), 'v', zeros(m.nd,nb), 'a', zeros(m.nd,nb));
    ib = 0;
end

cnt = struct('newton', 0, 'maxit', 0, 'cuts', 0, 'ev_sub', 0, 'steps', 0);
status = 'completed'; failed = false; t = 0; Lfail = NaN; mode = '';
for k = 1:nsteps
    t0 = (k-1)*dt;
    [q1, v1, a1, ok, nit, nc] = advance(m, q, v, a, t0, dt, W, o, 0);
    cnt.newton = cnt.newton + nit; cnt.cuts = cnt.cuts + nc;
    cnt.maxit = max(cnt.maxit, nit);
    if ~ok, status = 'no_convergence'; break; end
    [lamS, lamR, Mb] = lambda_of(m, q1, BSinv, capinv);
    L1 = crit(lamS, lamR, useS, useR);
    t1 = k*dt;
    if o.stop && L1 >= 1
        [q1, v1, a1, t1, ok, failed, Lfail, nsub, nitE, st] = localise(m, q, v, a, t0, dt, W, o, BSinv, capinv, useS, useR);
        cnt.ev_sub = cnt.ev_sub + nsub; cnt.newton = cnt.newton + nitE;
        if ~ok, status = st; break; end
        [lamS, lamR, Mb] = lambda_of(m, q1, BSinv, capinv);
        L1 = crit(lamS, lamR, useS, useR);
    end
    q = q1; v = v1; a = a1; t = t1; cnt.steps = k;

    % running statistics
    env = max(env, lamS);
    [ls, is] = max(lamS);
    if ls > lsmax, lsmax = ls; ism = is; tsm = t; end
    if lamR > lrmax, lrmax = lamR; trm = t; end
    if abs(Mb) > mbmax, mbmax = abs(Mb); end
    if isnan(tcross) && L1 >= 1
        tcross = t - dt*(L1 - 1)/max(L1 - Lprev, eps);
    end
    Lprev = L1;
    if t > o.stat_t0
        x = q(1:3:end); S.cnt = S.cnt + 1;
        S.sx = S.sx + x; S.sxx = S.sxx + x.^2; S.sy = S.sy + q(2:3:end);
        S.xmin = min(S.xmin, x); S.xmax = max(S.xmax, x);
        S.sl = S.sl + lamS; S.sll = S.sll + lamS.^2;
        S.sth = S.sth + q(3*n); S.sthth = S.sthth + q(3*n)^2;
        S.smb = S.smb + Mb; S.smbmb = S.smbmb + Mb^2;
    end
    if o.hist_every > 0 && (mod(k, o.hist_every) == 0 || failed)
        nh = nh + 1; wt = v35_wind_at(W, t);
        H(nh,:) = single([t, q(3*n-2), q(3*n-1), q(3*n), ls, is, Mb, wt(1)]);
    end
    if tailOn && (mod(k, o.tail_every) == 0 || failed)
        ib = mod(ib, nb) + 1;
        TB.t(ib) = t; TB.q(:,ib) = q; TB.v(:,ib) = v; TB.a(:,ib) = a;
    end
    if failed, break; end
end

out = struct();
out.status = status; out.ok = strcmp(status, 'completed');
out.failed = failed; out.t_end = t; out.t_fail = NaN; out.lam_fail = NaN;
out.node_fail = NaN; out.s_fail = NaN; out.mode = '';
if failed
    [lamS, lamR] = lambda_of(m, q, BSinv, capinv);
    [ls, is] = max(lamS);
    out.t_fail = t; out.lam_fail = Lfail;
    if useR && (~useS || lamR >= ls)
        out.mode = 'root'; out.node_fail = 1; out.s_fail = 0;
    else
        out.mode = 'stem'; out.node_fail = is; out.s_fail = m.s(is);
    end
end
out.q = q; out.v = v; out.a = a;
out.lam_stem_max = lsmax; out.node_stem_max = ism; out.s_stem_max = m.s(ism); out.t_stem_max = tsm;
out.lam_root_max = lrmax; out.t_root_max = trm; out.Mbase_absmax = mbmax;
out.t_first_passage = tcross;
out.lam_env = env;
c = max(S.cnt, 1);
out.stats = struct('n', S.cnt, 'x_mean', S.sx/c, 'x_std', sqrt(max(S.sxx/c - (S.sx/c).^2, 0)), ...
    'y_mean', S.sy/c, 'x_min', S.xmin, 'x_max', S.xmax, ...
    'lam_mean', S.sl/c, 'lam_rms', sqrt(S.sll/c), ...
    'th_tip_mean', S.sth/c, 'th_tip_std', sqrt(max(S.sthth/c - (S.sth/c)^2, 0)), ...
    'Mbase_mean', S.smb/c, 'Mbase_std', sqrt(max(S.smbmb/c - (S.smb/c)^2, 0)));
out.hist = H(1:nh,:);
if tailOn
    ord = mod(ib + (1:nb) - 1, nb) + 1; ord = ord(~isnan(TB.t(ord)));
    out.tail = struct('t', TB.t(ord), 'q', TB.q(:,ord), 'v', TB.v(:,ord), 'a', TB.a(:,ord));
end
out.counters = cnt; out.dt = dt;
end

%% ------------------------------------------------------------------------
function o = defaults(o, W, m)
d = struct('T', W.T, 'stop', true, 'tol', 1e-9, 'max_iter', 12, 'max_cut', 5, ...
    'ev_tol', 2e-3, 'ev_depth', 14, 'ev_guard', 1.05, 'stat_t0', W.T_ramp, ...
    'hist_every', 0, 'tail_T', 0, 'tail_every', 1);
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(o, fn{i}) || isempty(o.(fn{i})), o.(fn{i}) = d.(fn{i}); end
end
if ~isfield(o, 'dt'), error('V35:newmark','opts.dt is required.'); end
o.T = min(o.T, W.T);
o.tol = max(o.tol, m.R_floor);
end

function L = crit(lamS, lamR, useS, useR)
L = -inf;
if useS, L = max(lamS); end
if useR, L = max(L, lamR); end
end

function [lamS, lamR, Mb] = lambda_of(m, q, BSinv, capinv)
th = q(3:3:end);
Me = m.EI.*diff(th)./m.h;
if m.P.clamped, Mb = Me(1); else, Mb = m.k_theta*(th(1) - m.theta_ref); end
lamS = abs([Mb; 0.5*(Me(1:end-1) + Me(2:end)); Me(end)]).*BSinv;
lamR = abs(Mb)*capinv;
end

function [q1, v1, a1, ok, nit, ncut] = advance(m, q, v, a, t, h, W, o, depth)
[q1, v1, a1, ok, nit] = nstep(m, q, v, a, t, h, W, o);
ncut = 0;
if ok || depth >= o.max_cut, return; end
[qm, vm, am, ok1, n1, c1] = advance(m, q, v, a, t, h/2, W, o, depth+1);
if ~ok1, ok = false; nit = nit + n1; ncut = c1 + 1; return; end
[q1, v1, a1, ok, n2, c2] = advance(m, qm, vm, am, t + h/2, h/2, W, o, depth+1);
nit = nit + n1 + n2; ncut = c1 + c2 + 1;
end

function [q, v, a, ok, it] = nstep(m, qn, vn, an, tn, h, W, o)
% One Newmark step: Newton with the consistent tangent; a step that raises
% the residual is halved (backtracking).
f = m.free;
a0c = 4/h^2; a1c = 2/h;
qp = qn + h*vn + (h^2/4)*an; vp = vn + (h/2)*an;
q = qn + h*vn;                              % predictor a_{n+1} = -a_n (stiff modes alternate)
wt = v35_wind_at(W, tn + h);
ok = false; nrOld = inf; dq = []; al = 1; qOld = q(f);
for it = 1:o.max_iter
    a = a0c*(q - qp); v = vp + (h/2)*a;
    [R, K] = v35_residual(m, q, v, a, wt, a0c, a1c);
    r = R(f); nr = norm(r, inf);
    if nr <= o.tol, ok = true; break; end
    if ~isfinite(nr), break; end
    if ~isempty(dq) && nr > nrOld
        if al <= 1/64, break; end
        al = al/2; q(f) = qOld + al*dq;           % backtrack
        continue;
    end
    if ~isempty(dq) && norm(al*dq, inf) < 1e-14*m.L && nr < 1e3*o.tol
        ok = true; break;                        % round-off floor
    end
    qOld = q(f); nrOld = nr; al = 1;
    dq = -v35_solve(K, r, m.damp, a1c);
    if any(~isfinite(dq)), break; end
    q(f) = qOld + dq;
end
a = a0c*(q - qp); v = vp + (h/2)*a;
end

function [q, v, a, t, ok, failed, L, nsub, nit, st] = localise(m, q, v, a, t0, H, W, o, BSinv, capinv, useS, useR)
% First passage inside [t0, t0+H]: integrate with sub-steps from the last
% accepted state, halving the sub-step whenever it crosses lambda = 1 by
% more than ev_tol. The wind record is unchanged.
te = t0 + H; t = t0; hs = H/2; hmin = H/2^o.ev_depth;
ok = true; failed = false; L = NaN; nsub = 0; nit = 0; st = 'completed';
while t < te - 1e-12*H
    hs = min(hs, te - t);
    [q1, v1, a1, okS, n1] = advance(m, q, v, a, t, hs, W, o, 0);
    nsub = nsub + 1; nit = nit + n1;
    if ~okS, ok = false; st = 'event_no_convergence'; return; end
    [lamS, lamR] = lambda_of(m, q1, BSinv, capinv);
    L1 = crit(lamS, lamR, useS, useR);
    if L1 < 1
        q = q1; v = v1; a = a1; t = t + hs;
    elseif L1 <= 1 + o.ev_tol
        q = q1; v = v1; a = a1; t = t + hs; failed = true; L = L1; return;
    elseif hs <= hmin*(1 + 1e-9)
        if L1 <= o.ev_guard
            q = q1; v = v1; a = a1; t = t + hs; failed = true; L = L1; return;
        end
        ok = false; st = 'event_guard_exceeded'; L = L1; return;
    else
        hs = hs/2;
    end
end
end
