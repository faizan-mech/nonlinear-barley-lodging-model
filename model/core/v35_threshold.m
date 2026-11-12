function R = v35_threshold(caseIds, set)
%V35_THRESHOLD Per-seed failure-speed search for a set of cases, in parallel.
%
% Bracket each survival-failure transition using geometric quarter points
% until its relative width is below set.tol. Completed records are saved
% for restart. The search assumes monotonicity and reports observed reversals.
%
% set fields (defaults in brackets):
%   N [61], dt [2.5e-3], seeds [1:24], tol [0.005], outDir, cacheDir,
%   Ncal [481], U50_ref [] (pilot median for the reference, sets the guess),
%   g [3.75] (guess U = Uqs0/(1+g*Iu) when U50_ref is empty), Umin [0.5],
%   Umax [60], max_hours [70], ropts (struct passed to V35_RECORD),
%   ropts_ref (extra options for case 'ref'), verbose [true]
%
% R(c) holds the per-seed thresholds and U50 statistics; a summary.mat is
% written for each completed case.

set = defaults(set);
C = v35_find_case(caseIds{1});
for i = 2:numel(caseIds), C(i) = v35_find_case(caseIds{i}); end
nc = numel(C);
t0 = tic;
logmsg(set, 'V35 threshold: %d cases, N = %d, dt = %.4g s, %d seeds', nc, set.N, set.dt, numel(set.seeds));

%% models and guesses
models = cell(nc,1); Ps = cell(nc,1); guess = zeros(nc,1); ropts = cell(nc,1);
gref = [];
for c = 1:nc
    cdir = fullfile(set.outDir, C(c).id);
    if ~exist(cdir, 'dir'), mkdir(cdir); end
    [models{c}, Ps{c}, S] = v35_case_model(C(c), set.N, set.cacheDir, set.Ncal);
    ro = set.ropts; ro.dt = set.dt;
    if strcmp(C(c).id, 'ref'), ro = merge(ro, set.ropts_ref); end
    ropts{c} = ro;
    mf = fullfile(cdir, 'meta.mat');
    if exist(mf, 'file')
        L = load(mf); guess(c) = L.meta.guess;
    else
        if isempty(gref), gref = guess_factor(set); end
        Uq0 = v35_qs_threshold(models{c}, 0);
        guess(c) = Uq0/(1 + gref*Ps{c}.wind.Iu);
        meta = struct('casedef', C(c), 'guess', guess(c), 'Uqs0', Uq0, 'g', gref, ...
            'N', set.N, 'dt', set.dt, 'seeds', set.seeds, 'tol', set.tol, ...
            'supports', S, 'f_grav', models{c}.f_grav, 'E', models{c}.E, ...
            'root_cap', models{c}.fail.root_cap);
        v35_save_atomic(mf, struct('meta', meta));
    end
    logmsg(set, '  %-14s guess %.3f m/s, f1 %.4f Hz, k_theta %.4f, k_M %.4f, k_W %.4f', ...
        C(c).id, guess(c), models{c}.f_grav(1), models{c}.k_theta, models{c}.k_M, models{c}.k_W);
end

%% searches (case x seed), replayed from saved records
ns = numel(set.seeds);
srch = repmat(struct('c', 0, 'seed', 0, 'lo', NaN, 'hi', NaN, 'ev', zeros(0,3), ...
    'rU', [], 'rres', [], 'rsub', [], 'done', false, 'status', 'active', ...
    'Ustar', NaN, 'nmono', 0), nc*ns, 1);
queue = zeros(0,2);                  % [search index, U]
k = 0;
for c = 1:nc
    for j = 1:ns
        k = k + 1;
        srch(k).c = c; srch(k).seed = set.seeds(j);
        [srch(k), newU] = advance(srch(k), guess(c), set, C(c).id);
        queue = [queue; [k*ones(numel(newU),1), newU(:)]]; %#ok<AGROW>
    end
end
logmsg(set, '  %d searches, %d already complete, %d records queued', numel(srch), sum([srch.done]), size(queue,1));

%% run
isOct = exist('OCTAVE_VERSION', 'builtin') ~= 0;
pool = [];
if ~isOct, pool = gcp('nocreate'); end
nDone = 0; lastLog = toc(t0); stopSubmit = false;
if isempty(pool)
    while ~isempty(queue)
        if toc(t0)/3600 > set.max_hours, break; end
        kk = queue(1,1); U = queue(1,2); queue(1,:) = [];
        c = srch(kk).c;
        res = v35_task(models{c}, Ps{c}, U, srch(kk).seed, ropts{c}, recfile(set, C(c).id, srch(kk).seed, U));
        [srch(kk), newU] = receive(srch(kk), U, res, guess(c), set, C(c).id);
        queue = [queue; [kk*ones(numel(newU),1), newU(:)]]; %#ok<AGROW>
        nDone = nDone + 1;
        if toc(t0) - lastLog > set.log_every, progress(set, srch, C, nDone, size(queue,1), 0, t0); lastLog = toc(t0); end
    end
else
    maxFly = 2*pool.NumWorkers;
    F = parallel.FevalFuture.empty; meta = zeros(0,2);
    while true
        if ~stopSubmit && toc(t0)/3600 > set.max_hours
            stopSubmit = true; logmsg(set, '  time limit reached: finishing running records');
        end
        while ~stopSubmit && numel(F) < maxFly && ~isempty(queue)
            kk = queue(1,1); U = queue(1,2); queue(1,:) = [];
            c = srch(kk).c;
            F(end+1) = parfeval(pool, @v35_task, 1, models{c}, Ps{c}, U, srch(kk).seed, ...
                ropts{c}, recfile(set, C(c).id, srch(kk).seed, U)); %#ok<AGROW>
            meta(end+1,:) = [kk, U]; %#ok<AGROW>
        end
        if isempty(F), break; end
        [i, res] = fetchNext(F);
        kk = meta(i,1); U = meta(i,2);
        F(i) = []; meta(i,:) = [];
        if ~isempty(res.message), logmsg(set, '  record error (%s s%d U %.4f): %s', C(srch(kk).c).id, srch(kk).seed, U, res.message); end
        [srch(kk), newU] = receive(srch(kk), U, res, guess(srch(kk).c), set, C(srch(kk).c).id);
        queue = [queue; [kk*ones(numel(newU),1), newU(:)]]; %#ok<AGROW>
        nDone = nDone + 1;
        if toc(t0) - lastLog > set.log_every, progress(set, srch, C, nDone, size(queue,1), numel(F), t0); lastLog = toc(t0); end
    end
end
progress(set, srch, C, nDone, size(queue,1), 0, t0);

%% summaries
R = struct('id', {}, 'complete', {}, 'Ustar', {}, 'stats', {});
for c = 1:nc
    ks = find([srch.c] == c);
    R(c).id = C(c).id;
    R(c).complete = all([srch(ks).done]);
    R(c).Ustar = [srch(ks).Ustar];
    if R(c).complete
        R(c).stats = summarise(set, C(c), srch(ks), models{c}, guess(c));
        logmsg(set, '  %-14s U50 = %.3f m/s [%.3f, %.3f], lognormal %.3f, non-monotone %d, excluded %d', ...
            C(c).id, R(c).stats.U50, R(c).stats.U50_ci, R(c).stats.U50_lognormal, ...
            R(c).stats.n_nonmonotone, R(c).stats.n_excluded);
    end
end
end

%% ------------------------------------------------------------------------
function set = defaults(set)
d = struct('N', 61, 'dt', 2.5e-3, 'seeds', 1:24, 'tol', 0.005, 'outDir', 'results/threshold', ...
    'cacheDir', 'cache', 'Ncal', 481, 'U50_ref', [], 'g', 3.75, 'Umin', 0.5, 'Umax', 60, ...
    'max_hours', 70, 'ropts', struct(), 'ropts_ref', struct(), 'verbose', true, ...
    'log_every', 300, 'logfile', '');
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(set, fn{i}), set.(fn{i}) = d.(fn{i}); end
end
end

function g = guess_factor(set)
% g in U50 (1 + g Iu) = Uqs0, from the pilot reference median when known.
if isempty(set.U50_ref), g = set.g; return; end
[mr, Pr] = v35_case_model(v35_find_case('ref'), set.N, set.cacheDir, set.Ncal);
Uq0 = v35_qs_threshold(mr, 0);
g = (Uq0/set.U50_ref - 1)/Pr.wind.Iu;
end

function f = recfile(set, id, seed, U)
f = fullfile(set.outDir, id, sprintf('s%03d_U%010.5f.mat', seed, U));
end

function [s, newU] = advance(s, g, set, id)
% Replay saved rounds; return the speeds of the next round still to run.
newU = [];
while ~s.done
    if isempty(s.rU)
        s.rU = propose(s, g, set); s.rres = nan(size(s.rU)); s.rsub = false(size(s.rU));
    end
    for j = 1:numel(s.rU)
        if isnan(s.rres(j))
            f = recfile(set, id, s.seed, s.rU(j));
            if exist(f, 'file')
                L = load(f); s.rres(j) = double(L.r.failed) + 2*double(~L.r.ok);
            end
        end
    end
    if any(isnan(s.rres))
        j = isnan(s.rres) & ~s.rsub;
        newU = s.rU(j); s.rsub(j) = true;
        return;
    end
    s = close_round(s, set);
end
end

function [s, newU] = receive(s, U, res, g, set, id)
j = find(abs(s.rU - U) < 1e-9, 1);
if isempty(j), newU = []; return; end
if strcmp(res.status, 'error') || ~res.ok
    s.rres(j) = 2;
else
    s.rres(j) = double(res.failed);
end
newU = [];
if all(~isnan(s.rres))
    s = close_round(s, set);
    [s, newU] = advance(s, g, set, id);
end
end

function s = close_round(s, set)
s.ev = [s.ev; [s.rU(:), double(s.rres(:) == 1), double(s.rres(:) ~= 2)]];
s.rU = []; s.rres = []; s.rsub = [];
if any(s.ev(:,3) == 0)
    s.done = true; s.status = 'numerical'; s.Ustar = NaN; return;
end
fl = s.ev(s.ev(:,2) == 1, 1); sv = s.ev(s.ev(:,2) == 0, 1);
if isempty(fl), s.hi = NaN; else, s.hi = min(fl); end
if isnan(s.hi)
    if isempty(sv), s.lo = NaN; else, s.lo = max(sv); end
    s.nmono = 0;
else
    b = sv(sv < s.hi);
    if isempty(b), s.lo = NaN; else, s.lo = max(b); end
    s.nmono = sum(sv > s.hi);
end
if ~isnan(s.lo) && ~isnan(s.hi) && s.hi/s.lo - 1 <= set.tol
    s.done = true; s.status = 'bracketed'; s.Ustar = sqrt(s.lo*s.hi);
elseif isnan(s.hi) && ~isnan(s.lo) && s.lo >= set.Umax
    s.done = true; s.status = 'censored_high'; s.Ustar = Inf;
elseif isnan(s.lo) && ~isnan(s.hi) && s.hi <= set.Umin
    s.done = true; s.status = 'censored_low'; s.Ustar = s.hi;
end
end

function U = propose(s, g, set)
r = 1.25;
if isnan(s.lo) && isnan(s.hi)
    U = g*[1/1.15, 1, 1.15];
elseif isnan(s.lo)
    U = s.hi*r.^(-3:-1);
elseif isnan(s.hi)
    U = s.lo*r.^(1:3);
else
    U = s.lo*(s.hi/s.lo).^((1:3)/4);
end
U = min(max(U, set.Umin), set.Umax);
U = unique(round(U*1e5)/1e5);
if ~isempty(s.ev), U = U(~ismember(U, s.ev(:,1))); end
if isempty(U)                        % bracket limits reached
    U = s.lo*(s.hi/s.lo)^0.5; U = round(U*1e5)/1e5;
end
end

function st = summarise(set, cs, S, m, g)
cdir = fullfile(set.outDir, cs.id);
Ustar = [S.Ustar];
ok = ~strcmp({S.status}, 'numerical');
st = v35_u50_stats(Ustar(ok));
st.Ustar = Ustar; st.seeds = [S.seed]; st.status = {S.status};
st.lo = [S.lo]; st.hi = [S.hi];
st.n_nonmonotone = sum([S.nmono] > 0); st.n_excluded = sum(~ok);
st.guess = g; st.casedef = cs; st.N = set.N; st.dt = set.dt; st.tol = set.tol;
% failure site and time of the lowest failing record of each seed
ns = numel(S); st.s_fail = nan(1,ns); st.t_fail = nan(1,ns); st.mode = cell(1,ns);
st.budget_share = nan(ns, 8);
for j = 1:ns
    if ~isnan(S(j).hi)
        L = load(recfile(set, cs.id, S(j).seed, S(j).hi)); r = L.r;
        st.s_fail(j) = r.s_fail; st.t_fail(j) = r.t_fail; st.mode{j} = r.mode;
        if ~isempty(r.budget), st.budget_share(j,:) = r.budget.share; end
    end
end
% digest of every record of the case
F = dir(fullfile(cdir, 's*_U*.mat'));
D = zeros(numel(F), 10);
for i = 1:numel(F)
    L = load(fullfile(cdir, F(i).name)); r = L.r;
    D(i,:) = [r.seed, r.U, r.failed, r.ok, r.t_fail, r.s_fail, r.lam_stem_max, ...
              r.s_stem_max, r.lam_root_max, r.wall];
end
st.records = D;
st.records_cols = {'seed','U','failed','ok','t_fail','s_fail','lam_stem_max','s_stem_max','lam_root_max','wall'};
st.f_grav = m.f_grav; st.k_theta = m.k_theta; st.k_M = m.k_M; st.k_W = m.k_W;
st.L = m.L; st.s = m.s;
v35_save_atomic(fullfile(cdir, 'summary.mat'), struct('st', st));
end

function progress(set, S, C, nDone, nQ, nRun, t0)
act = ~[S.done];
cs = unique([S(act).c]);
names = '';
for c = cs(:).', names = [names ' ' C(c).id]; end %#ok<AGROW>
logmsg(set, '  [%.2f h] records done %d, running %d, queued %d | searches done %d/%d | active cases:%s', ...
    toc(t0)/3600, nDone, nRun, nQ, sum(~act), numel(S), names);
end

function logmsg(set, varargin)
msg = sprintf(varargin{:});
if set.verbose, fprintf('%s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg); end
if ~isempty(set.logfile)
    fid = fopen(set.logfile, 'a'); fprintf(fid, '%s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg); fclose(fid);
end
end

function a = merge(a, b)
fn = fieldnames(b);
for i = 1:numel(fn), a.(fn{i}) = b.(fn{i}); end
end
