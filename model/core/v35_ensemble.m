function D = v35_ensemble(jobs, set)
%V35_ENSEMBLE Run a fixed list of records in parallel, resumable.
%
% jobs: struct array with fields cid (case id), N, dt, U, seed, file, and
% optional ropts (struct for V35_RECORD, e.g. stop = false, hist_every).
% set: cacheDir, Ncal, max_hours, verbose, logfile.
% Records already saved are skipped. Returns a digest per job.
if ~isfield(set, 'max_hours'), set.max_hours = 70; end
if ~isfield(set, 'verbose'), set.verbose = true; end
if ~isfield(set, 'logfile'), set.logfile = ''; end
if ~isfield(set, 'Ncal'), set.Ncal = 481; end
t0 = tic;
% models by (case, N)
keys = arrayfun(@(j) sprintf('%s|%d', j.cid, j.N), jobs, 'UniformOutput', false);
[uk, ~, ik] = unique(keys);
models = cell(numel(uk),1); Ps = cell(numel(uk),1);
for i = 1:numel(uk)
    j = find(ik == i, 1);
    [models{i}, Ps{i}] = v35_case_model(jobs(j).cid, jobs(j).N, set.cacheDir, set.Ncal);
end
todo = find(arrayfun(@(j) ~exist(j.file, 'file'), jobs));
logmsg(set, 'V35 ensemble: %d records, %d to run, %d models', numel(jobs), numel(todo), numel(uk));
% expensive records first (fine mesh, small step)
cost = arrayfun(@(j) j.N/j.dt, jobs(todo));
[~, o] = sort(cost, 'descend'); todo = todo(o);
isOct = exist('OCTAVE_VERSION', 'builtin') ~= 0;
pool = [];
if ~isOct, pool = gcp('nocreate'); end
ropt = @(j) setdt(j);
nDone = 0; lastLog = 0;
if isempty(pool)
    for q = todo(:).'
        if toc(t0)/3600 > set.max_hours, break; end
        j = jobs(q);
        res = v35_task(models{ik(q)}, Ps{ik(q)}, j.U, j.seed, ropt(j), j.file);
        nDone = nDone + 1;
        if ~isempty(res.message), logmsg(set, '  error %s: %s', j.file, res.message); end
        if toc(t0) - lastLog > 300, logmsg(set, '  [%.2f h] %d/%d', toc(t0)/3600, nDone, numel(todo)); lastLog = toc(t0); end
    end
else
    F = parallel.FevalFuture.empty; meta = zeros(0,1); next = 1; maxFly = 2*pool.NumWorkers;
    while true
        while next <= numel(todo) && numel(F) < maxFly && toc(t0)/3600 <= set.max_hours
            q = todo(next); next = next + 1; j = jobs(q);
            F(end+1) = parfeval(pool, @v35_task, 1, models{ik(q)}, Ps{ik(q)}, j.U, j.seed, ropt(j), j.file); %#ok<AGROW>
            meta(end+1) = q; %#ok<AGROW>
        end
        if isempty(F), break; end
        [i, res] = fetchNext(F);
        if ~isempty(res.message), logmsg(set, '  error %s: %s', jobs(meta(i)).file, res.message); end
        F(i) = []; meta(i) = [];
        nDone = nDone + 1;
        if toc(t0) - lastLog > 300, logmsg(set, '  [%.2f h] %d/%d, running %d', toc(t0)/3600, nDone, numel(todo), numel(F)); lastLog = toc(t0); end
    end
end
% digest
D = repmat(struct('file', '', 'done', false, 'failed', NaN, 'ok', false, 'lam_stem_max', NaN, ...
    's_stem_max', NaN, 't_first_passage', NaN, 't_fail', NaN, 's_fail', NaN, 'mode', ''), numel(jobs), 1);
for q = 1:numel(jobs)
    D(q).file = jobs(q).file;
    if exist(jobs(q).file, 'file')
        L = load(jobs(q).file); r = L.r;
        D(q).done = true; D(q).failed = r.failed; D(q).ok = r.ok; D(q).lam_stem_max = r.lam_stem_max;
        D(q).s_stem_max = r.s_stem_max; D(q).t_first_passage = r.t_first_passage;
        D(q).t_fail = r.t_fail; D(q).s_fail = r.s_fail; D(q).mode = r.mode;
    end
end
logmsg(set, 'V35 ensemble: %d/%d records complete (%.2f h)', sum([D.done]), numel(jobs), toc(t0)/3600);
end

function o = setdt(j)
if isfield(j, 'ropts') && ~isempty(j.ropts), o = j.ropts; else, o = struct(); end
o.dt = j.dt;
end

function logmsg(set, varargin)
msg = sprintf(varargin{:});
if set.verbose, fprintf('%s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg); end
if ~isempty(set.logfile)
    fid = fopen(set.logfile, 'a'); fprintf(fid, '%s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg); fclose(fid);
end
end
