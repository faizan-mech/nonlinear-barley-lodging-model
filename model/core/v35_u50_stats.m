function S = v35_u50_stats(Ustar, nboot)
%V35_U50_STATS Median failure speed from per-seed thresholds, with intervals.
%
% Returns the median, a deterministic percentile-bootstrap 95 % interval,
% a lognormal fit and its 10 % and 90 % failure speeds. The probability
% interpretation assumes a monotone failure transition for each seed.
if nargin < 2, nboot = 4000; end
U = Ustar(:); U = U(~isnan(U));
n = numel(U);
S = struct('n', n, 'n_censored', sum(isinf(U)), 'U50', NaN, 'U50_ci', [NaN NaN], ...
    'mu_log', NaN, 'sd_log', NaN, 'U50_lognormal', NaN, 'U50_lognormal_ci', [NaN NaN], ...
    'U10', NaN, 'U90', NaN, 'U10_emp', NaN, 'U90_emp', NaN);
if n == 0, return; end
S.U50 = median(U);
r = v35_uniform(20260930, 7, n*nboot);
idx = min(n, 1 + floor(n*reshape(r, n, nboot)));
med = sort(median(U(idx), 1));
S.U50_ci = [med(max(1, round(0.025*nboot))), med(round(0.975*nboot))];
Uf = U(isfinite(U));
if numel(Uf) == n && n >= 3
    lu = log(Uf); S.mu_log = mean(lu); S.sd_log = std(lu);
    S.U50_lognormal = exp(S.mu_log);
    se = S.sd_log/sqrt(n);
    S.U50_lognormal_ci = exp(S.mu_log + [-1 1]*1.96*se);
    S.U10 = exp(S.mu_log - 1.2816*S.sd_log); S.U90 = exp(S.mu_log + 1.2816*S.sd_log);
end
Us = sort(U);
S.U10_emp = quant(Us, 0.1); S.U90_emp = quant(Us, 0.9);
end

function x = quant(Us, p)
n = numel(Us); h = (n-1)*p + 1; lo = floor(h); hi = ceil(h);
if hi == lo, x = Us(lo); else, x = Us(lo) + (h - lo)*(Us(hi) - Us(lo)); end
end
