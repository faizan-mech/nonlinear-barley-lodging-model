function S = v35_supports(cs, P, cacheDir, Ncal)
%V35_SUPPORTS Support constants of case cs (calibrated once, then cached).
% The reference calibration is stored as calib_ref.mat; cases with their
% own calibration as calib_<sid>.mat. Files are written atomically, so
% concurrent jobs can share the cache.
if nargin < 4 || isempty(Ncal), Ncal = 481; end
Sref = load_or_make('ref', 'both', v35_params(), [], cacheDir, Ncal);
switch cs.calib
    case 'none'
        S = Sref;
        S.k_theta = Sref.k_theta*P.k_theta_factor;
        S.rule = 'none';
    otherwise
        S = load_or_make(cs.sid, cs.calib, P, Sref, cacheDir, Ncal);
end
end

function S = load_or_make(sid, rule, P, Sref, cacheDir, Ncal)
f = fullfile(cacheDir, ['calib_' sid '.mat']);
if exist(f, 'file')
    L = load(f); S = L.S;
    if S.Ncal == Ncal, return; end
end
S = v35_calibrate(P, rule, Sref, Ncal);
S.sid = sid;
v35_save_atomic(f, struct('S', S));
end
