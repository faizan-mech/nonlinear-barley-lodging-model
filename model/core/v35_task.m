function res = v35_task(m, P, U, seed, ropts, outfile)
%V35_TASK Worker entry: run one record, save it atomically, return a digest.
% Errors are caught and returned, so a failing record never stops a pool.
res = struct('U', U, 'seed', seed, 'ok', false, 'failed', false, ...
    'status', '', 'file', outfile, 'wall', NaN, 'message', '');
try
    if ~isempty(outfile) && exist(outfile, 'file')
        L = load(outfile); r = L.r;
    else
        r = v35_record(m, P, U, seed, ropts);
        if ~isempty(outfile), v35_save_atomic(outfile, struct('r', r)); end
    end
    res.ok = r.ok; res.failed = r.failed; res.status = r.status; res.wall = r.wall;
catch err
    res.status = 'error'; res.message = err.message;
end
end
