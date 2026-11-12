function [q, info] = v35_equilibrium(m, Ueq, q0, opts)
%V35_EQUILIBRIUM Mean-wind equilibrium by continuation in the load factor.
%
% The load factor p = (U/Ueq)^2 runs from 0 (gravity equilibrium) to 1 with
% adaptive steps and a secant predictor. The ear sees Ueq, the culm
% alpha(s)*Ueq; the vertical component has zero mean.
if nargin < 3 || isempty(q0), q0 = m.q_grav; end
if nargin < 4, opts = struct(); end
info = struct('ok', false, 'steps', 0, 'newton', 0, 'cuts', 0);
q = q0;
if Ueq <= 0, info.ok = true; return; end
p = 0; dp = 0.25; qprev = []; pprev = [];
while p < 1
    pn = min(1, p + dp);
    if isempty(qprev), qg = q; else, qg = q + (q - qprev)*(pn - p)/(p - pprev); end
    [qn, ok, it] = v35_static(m, qg, [sqrt(pn)*Ueq, 0], opts);
    info.newton = info.newton + it;
    if ~ok && ~isempty(qprev)          % retry without predictor
        [qn, ok, it] = v35_static(m, q, [sqrt(pn)*Ueq, 0], opts);
        info.newton = info.newton + it;
    end
    if ok
        qprev = q; pprev = p; q = qn; p = pn; info.steps = info.steps + 1;
        if it <= 5, dp = min(0.5, 1.5*dp); end
    else
        dp = dp/2; info.cuts = info.cuts + 1;
        if dp < 1e-5, return; end
    end
end
info.ok = true;
end
