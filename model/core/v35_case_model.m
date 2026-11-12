function [m, P, S] = v35_case_model(cs, N, cacheDir, Ncal)
%V35_CASE_MODEL Parameters, supports and prepared model of case cs on N nodes.
if nargin < 4, Ncal = []; end
if ischar(cs), cs = v35_find_case(cs); end
P = v35_apply_variant(v35_params(), cs.ov);
S = v35_supports(cs, P, cacheDir, Ncal);
P.k_theta = S.k_theta; P.k_M = S.k_M; P.k_W = S.k_W;
P.k_theta_factor = 1;                 % already applied in S
m = v35_build(P, N);
m = v35_prepare(m);
m.case_id = cs.id;
end
