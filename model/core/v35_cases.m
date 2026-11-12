function C = v35_cases()
%V35_CASES Reference, inflow, trait and support cases.
%
% Fields: id, group, desc, ov (V35_APPLY_VARIANT overrides), calib and sid.
% calib: 'none' keeps reference supports; 'kM' recalibrates canopy support;
% 'both' recalibrates root and canopy supports. sid is the calibration-cache key.

C = struct('id', {}, 'group', {}, 'desc', {}, 'ov', {}, 'calib', {}, 'sid', {});
C = add(C, 'ref', 'ref', 'Reference: Iu 0.40, Lu = L, von Karman, crossflow ear', {}, 'both');

% inflow
for x = [0.2 0.3 0.5 0.6]
    C = add(C, sprintf('iu%03d', round(100*x)), 'inflow', sprintf('Iu = %.2f', x), {'wind.Iu', x}, 'none');
end
for x = [0.25 0.5 1.5 2.0]
    C = add(C, sprintf('lu%03d', round(100*x)), 'inflow', sprintf('Lu = %.2f L', x), {'wind.LuL', x}, 'none');
end
C = add(C, 'kaimal', 'inflow', 'Kaimal spectrum', {'wind.spectrum', 'kaimal'}, 'none');
C = add(C, 'vertical', 'inflow', 'Correlated vertical fluctuation', {'wind.vertical', true}, 'none');

% ear load law
C = add(C, 'ear_upright', 'ear', 'Ear load at upright area and arm (constant load)', {'ear.law', 'upright'}, 'none');
C = add(C, 'ear_projection', 'ear', 'Ear load by projected area (bound)', {'ear.law', 'projection'}, 'none');

% ear mass pathways
C = add(C, 'ear_massless', 'ear', 'Aerodynamic ear only', ...
    {'ear.s_inertia', 0, 'ear.s_weight_force', 0, 'ear.s_weight_couple', 0}, 'none');
C = add(C, 'ear_inertia', 'ear', 'Ear inertia only', ...
    {'ear.s_weight_force', 0, 'ear.s_weight_couple', 0}, 'none');
C = add(C, 'ear_nocouple', 'ear', 'Inertia and weight, no weight couple', ...
    {'ear.s_weight_couple', 0}, 'none');
for f = [0.8 1.2]
    tag = sprintf('%03d', round(100*f));
    C = add(C, ['ear_si' tag], 'ear', sprintf('Ear inertia x%.1f', f), {'ear.s_inertia', f}, 'none');
    C = add(C, ['ear_swf' tag], 'ear', sprintf('Ear weight force x%.1f', f), {'ear.s_weight_force', f}, 'none');
    C = add(C, ['ear_swc' tag], 'ear', sprintf('Ear weight couple x%.1f', f), {'ear.s_weight_couple', f}, 'none');
    C = add(C, ['ear_area' tag], 'ear', sprintf('Ear areas x%.1f', f), {'ear.A_up*', f, 'ear.A_hor*', f}, 'none');
    C = add(C, ['ear_mass' tag], 'ear', sprintf('Ear mass and inertia x%.1f', f), {'ear.m*', f, 'ear.J_G*', f}, 'none');
end

% culm properties (+/-20 %)
for f = [0.8 1.2]
    tag = sprintf('%03d', round(100*f));
    C = add(C, ['bs' tag], 'culm', sprintf('Strength x%.1f', f), {'scale_BS', f}, 'none');
    C = add(C, ['do_t_bs' tag], 'culm', sprintf('Do x%.1f, t fixed, BS fixed', f), {'scale_Do', f}, 'none');
    C = add(C, ['do_t_z' tag], 'culm', sprintf('Do x%.1f, t fixed, BS ~ Z', f), {'scale_Do', f, 'BS_scales_with_Z', true}, 'none');
    C = add(C, ['do_tp_bs' tag], 'culm', sprintf('Do x%.1f, t ~ Do, BS fixed', f), {'scale_Do', f, 'scale_t', f}, 'none');
    C = add(C, ['do_tp_z' tag], 'culm', sprintf('Do x%.1f, t ~ Do, BS ~ Z', f), {'scale_Do', f, 'scale_t', f, 'BS_scales_with_Z', true}, 'none');
    C = add(C, ['t_out' tag], 'culm', sprintf('Wall x%.1f, outer radius fixed', f), {'scale_t', f}, 'none');
    C = add(C, ['t_in' tag], 'culm', sprintf('Wall x%.1f, inner radius fixed', f), {'scale_t', f, 'wall_keep', 'inner'}, 'none');
    C = add(C, ['e' tag], 'culm', sprintf('Modulus x%.1f', f), {'scale_E', f}, 'none');
end
C = add(C, 'zeta004', 'culm', 'Damping ratio 0.04', {'zeta', 0.04}, 'none');
C = add(C, 'zeta012', 'culm', 'Damping ratio 0.12', {'zeta', 0.12}, 'none');

% alternative constructions and supports
P0 = v35_params(); L = P0.L;
C = add(C, 'bs_anchors', 'model', 'Strength by point anchors 0.10/0.40/0.70 L', ...
    {'BS_s', [0 0.10 0.40 0.70 1]*L, 'BS_v', [0.075 0.075 0.042 0.021 0.021]}, 'none');
C = add(C, 'ei_strength', 'model', 'Lower-culm EI proportional to strength (upper bound)', ...
    {'EI_follows_strength', true}, 'both');
C = add(C, 'winkler', 'model', 'Lateral (Winkler) canopy support', {'support', 'winkler'}, 'kM');
C = add(C, 'ktheta_N', 'model', 'Root spring k_theta/N (N = 4.4)', {'k_theta_factor', 1/4.4}, 'kM');
C = add(C, 'win_up', 'model', 'Window shifted +0.05 L', {'win_lo', 0.35, 'win_hi', 0.675}, 'kM');
C = add(C, 'win_dn', 'model', 'Window shifted -0.05 L', {'win_lo', 0.25, 'win_hi', 0.575}, 'kM');
C = add(C, 'win_d015', 'model', 'Window edge width 0.015 L', {'win_d', 0.015}, 'kM');
C = add(C, 'win_d060', 'model', 'Window edge width 0.06 L', {'win_d', 0.06}, 'kM');

% root anchorage (capacity shared by N shoots; stem criterion off)
for tau = [1 2 3 4 5]
    C = add(C, sprintf('root_tau%d', tau), 'root', sprintf('Root failure, tau_s = %d kPa', tau), ...
        {'fail.stem', false, 'fail.root_share', true, 'fail.tau_s', 1e3*tau}, 'none');
end
C = add(C, 'root_kth080', 'root', 'Root failure, tau_s 4 kPa, k_theta x0.8', ...
    {'fail.stem', false, 'fail.root_share', true, 'k_theta_factor', 0.8}, 'none');
C = add(C, 'root_kth120', 'root', 'Root failure, tau_s 4 kPa, k_theta x1.2', ...
    {'fail.stem', false, 'fail.root_share', true, 'k_theta_factor', 1.2}, 'none');
C = add(C, 'root_stem', 'root', 'First of root (M_R/N, 4 kPa) or stem failure', ...
    {'fail.root_share', true}, 'none');

% soil strength +/-20 % about 4 kPa
for tau = [3.2 4.8]
    C = add(C, sprintf('root_tau%dp%d', floor(tau), round(10*(tau - floor(tau)))), 'root', ...
        sprintf('Root failure, tau_s = %.1f kPa', tau), ...
        {'fail.stem', false, 'fail.root_share', true, 'fail.tau_s', 1e3*tau}, 'none');
end
% ear: projected area with the drag arm held upright (isolates area from arm)
C = add(C, 'ear_area_only', 'ear', 'Ear load by projected area, arm at its upright value', {'ear.law', 'area'}, 'none');
% strength profiles through the 20 cm section means of Berry et al. (2006),
% anchored at the section midpoints 0.10, 0.415 and 0.73 m
sm = [0.10 0.415 0.73]; bm = [0.075 0.042 0.021];
C = add(C, 'bs_mid', 'model', 'Strength by point anchors at the 20 cm section midpoints (linear)', ...
    {'BS_s', [0 sm L], 'BS_v', [bm(1) bm bm(3)]}, 'none');
sg = unique([linspace(0, L, 167), sm]);
bp = pchip(sm, bm, min(max(sg, sm(1)), sm(end)));
C = add(C, 'bs_pchip', 'model', 'Strength by monotone cubic (PCHIP) through the section midpoints', ...
    {'BS_s', sg, 'BS_v', bp}, 'none');
cf = polyfit(sm, log(bm), 1);                  % BS = B0 exp(-s/l), l = 0.495 m
C = add(C, 'bs_exp', 'model', 'Strength by exponential fit to the section midpoints', ...
    {'BS_s', sg, 'BS_v', exp(polyval(cf, sg))}, 'none');

ids = {C.id};
if numel(unique(ids)) ~= numel(ids), error('V35:cases','Duplicate case ids.'); end
end

function C = add(C, id, group, desc, ov, calib)
k = numel(C) + 1;
C(k).id = id; C(k).group = group; C(k).desc = desc; C(k).ov = ov; C(k).calib = calib;
if strcmp(calib, 'none'), C(k).sid = 'ref'; else, C(k).sid = id; end
end
