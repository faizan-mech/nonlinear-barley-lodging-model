function r = v35_postures(m, P, U, seed, opts)
%V35_POSTURES Swept region of the centreline over one record (no stopping).
%
% Continue the elastic trajectory after capacity crossing. Sample centrelines
% after the ramp and return the minimum and maximum x at each y = opts.ygrid*L.
%
% opts: dt (solver step), every (snapshot interval, s; default 0.02),
% ygrid (levels as fractions of L; default 0:0.005:1).

if ~isfield(opts, 'every'), opts.every = 0.02; end
if ~isfield(opts, 'ygrid'), opts.ygrid = 0:0.005:1; end
tw = tic; L = m.L; n = m.n;
W = v35_wind(P, U, seed);
r = struct('U', U, 'seed', seed, 'n', n, 'dt', opts.dt, 'ok', false, 'status', '', ...
    'y_L', opts.ygrid(:), 'xlo_L', [], 'xhi_L', [], 'x_mean_L', [], 'y_mean_L', [], ...
    'tip_L', [], 't', [], 'post_min_L', [], 'post_max_L', [], 'reach_max', NaN, ...
    'lam_stem_max', NaN, 't_first_passage', NaN, 'wall', NaN);
[qe, ie] = v35_equilibrium(m, W.Ueq);
if ~ie.ok, r.status = 'equilibrium_failed'; r.wall = toc(tw); return; end
out = [];
for fct = [1 0.5 0.25]
    dt = opts.dt*fct;
    o = struct('dt', dt, 'stop', false, 'tail_T', W.T + 1, 'tail_every', max(1, round(opts.every/dt)));
    out = v35_newmark(m, qe, W, o);
    if out.ok, break; end
end
r.status = out.status; r.ok = out.ok;
if ~out.ok, r.wall = toc(tw); return; end
k = out.tail.t > W.T_ramp;
t = out.tail.t(k); Q = out.tail.q(:, k);
X = Q(1:3:end, :)/L; Y = Q(2:3:end, :)/L;          % n x ns, normalised
yg = opts.ygrid(:); ny = numel(yg); ns = numel(t);
xlo = inf(ny,1); xhi = -inf(ny,1);
for j = 1:ns
    xa = X(1:end-1,j).'; xb = X(2:end,j).'; ya = Y(1:end-1,j).'; yb = Y(2:end,j).';
    % levels crossed by each segment (ny x n-1)
    lo = min(ya, yb); hi = max(ya, yb);
    hit = yg >= lo & yg <= hi;
    dy = yb - ya; dy(dy == 0) = eps;
    xs = xa + (yg - ya)./dy.*(xb - xa);
    xs(~hit) = NaN;
    xlo = min(xlo, min(xs, [], 2)); xhi = max(xhi, max(xs, [], 2));
end
xlo(~isfinite(xlo)) = NaN; xhi(~isfinite(xhi)) = NaN;
[~, jmin] = min(X(end,:)); [~, jmax] = max(X(end,:));
sL = m.s(:)/L;
reach = sqrt(X.^2 + Y.^2)./max(sL, eps); reach(1,:) = 0;
r.xlo_L = xlo; r.xhi_L = xhi;
r.x_mean_L = mean(X, 2); r.y_mean_L = mean(Y, 2);
r.tip_L = single([X(end,:).', Y(end,:).']); r.t = single(t(:));
r.post_min_L = [X(:,jmin), Y(:,jmin)]; r.post_max_L = [X(:,jmax), Y(:,jmax)];
r.reach_max = max(reach(:));
r.lam_stem_max = out.lam_stem_max; r.t_first_passage = out.t_first_passage;
r.wall = toc(tw);
end
