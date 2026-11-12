function r = v35_record(m, P, U, seed, opts)
%V35_RECORD One stochastic record: mean-load equilibrium, 2 s ramp, transient.
%
% r = v35_record(m, P, U, seed, opts), with m from V35_PREPARE and P the
% case parameters (inflow). opts: dt (solver step), stop (first passage,
% default true), retries (dt factors, default [1 0.5 0.25]), hist_every,
% tail_T, tail_every, budget (default true). Returns a compact summary.

if nargin < 5, opts = struct(); end
if ~isfield(opts, 'retries'), opts.retries = [1 0.5 0.25]; end
if ~isfield(opts, 'budget'), opts.budget = true; end
tw = tic;
W = v35_wind(P, U, seed);
r = struct('U', U, 'seed', seed, 'Ueq', W.Ueq, 'Iu', W.Iu, 'Lu', W.Lu, ...
    'spectrum', W.spectrum, 'n', m.n, 'dt', opts.dt, 'dt_used', NaN, ...
    'status', '', 'ok', false, 'failed', false, 'mode', '', 't_fail', NaN, ...
    'lam_fail', NaN, 'node_fail', NaN, 's_fail', NaN, 'lam_stem_max', NaN, ...
    's_stem_max', NaN, 't_stem_max', NaN, 'lam_root_max', NaN, 't_first_passage', NaN, ...
    'Mbase_absmax', NaN, 'lam_eq', NaN, 's_eq', NaN, 'tip_eq', [NaN NaN], ...
    'lam_env', [], 'stats', [], 'hist', [], 'budget', [], 'tail_budget', [], ...
    'counters', [], 'retries', 0, 'r_uw', W.r_uw_realised, 'wall', NaN);

[qe, ie] = v35_equilibrium(m, W.Ueq);
if ~ie.ok
    r.status = 'equilibrium_failed'; r.wall = toc(tw); return;
end
[Mn, Mb] = v35_moments(m, qe);
lamE = abs(Mn)./m.BS; [le, ie2] = max(lamE);
r.lam_eq = le; r.s_eq = m.s(ie2); r.tip_eq = [qe(3*m.n-2), qe(3*m.n-1)];
if isfield(opts, 'stop') && ~opts.stop, stopping = false; else, stopping = true; end
Lroot = abs(Mb)/m.fail.root_cap;
if stopping && ((m.fail.stem && le >= 1) || (m.fail.root && Lroot >= 1))
    r.status = 'completed'; r.ok = true; r.failed = true; r.t_fail = 0;
    if m.fail.root && Lroot >= max(le*m.fail.stem, 1)
        r.mode = 'root'; r.node_fail = 1; r.s_fail = 0; r.lam_fail = Lroot;
    else
        r.mode = 'stem'; r.node_fail = ie2; r.s_fail = m.s(ie2); r.lam_fail = le;
    end
    r.lam_stem_max = le; r.s_stem_max = m.s(ie2); r.lam_root_max = Lroot;
    r.wall = toc(tw); return;
end

o = opts; out = [];
for j = 1:numel(opts.retries)
    o.dt = opts.dt*opts.retries(j);
    if isfield(opts, 'hist_every') && opts.hist_every > 0
        o.hist_every = max(1, round(opts.hist_every/opts.retries(j)));
    end
    if isfield(opts, 'tail_every') && opts.tail_every > 0
        o.tail_every = max(1, round(opts.tail_every/opts.retries(j)));
    end
    out = v35_newmark(m, qe, W, o);
    r.retries = j - 1;
    if out.ok, break; end
end
r.dt_used = o.dt; r.status = out.status; r.ok = out.ok;
r.failed = out.failed; r.mode = out.mode; r.t_fail = out.t_fail; r.lam_fail = out.lam_fail;
r.node_fail = out.node_fail; r.s_fail = out.s_fail;
r.lam_stem_max = out.lam_stem_max; r.s_stem_max = out.s_stem_max; r.t_stem_max = out.t_stem_max;
r.lam_root_max = out.lam_root_max; r.t_first_passage = out.t_first_passage;
r.Mbase_absmax = out.Mbase_absmax;
r.lam_env = single(out.lam_env);
st = out.stats; fn = fieldnames(st);
for i = 1:numel(fn), if ~strcmp(fn{i}, 'n'), st.(fn{i}) = single(st.(fn{i})); end, end
r.stats = st;
r.hist = out.hist;
r.counters = out.counters;
if out.ok && out.failed && opts.budget
    k = out.node_fail;
    r.budget = v35_budget(m, out.q, out.v, out.a, v35_wind_at(W, out.t_end), k);
    if isfield(out, 'tail')
        T = out.tail; nb = numel(T.t); Mi = zeros(nb, 8); tot = zeros(nb,1); Mnk = zeros(nb,1);
        for i = 1:nb
            B = v35_budget(m, T.q(:,i), T.v(:,i), T.a(:,i), v35_wind_at(W, T.t(i)), k);
            Mi(i,:) = B.M; tot(i) = B.total; Mnk(i) = B.M_node;
        end
        r.tail_budget = struct('t', T.t(:), 'M', Mi, 'total', tot, 'M_node', Mnk, ...
            'BS', m.BS(k), 'names', {B.names});
    end
end
r.wall = toc(tw);
end
