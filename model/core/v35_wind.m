function W = v35_wind(P, U, seed)
%V35_WIND Synthesise one approximately Gaussian inflow record at ear height.
%
% u(t) = U + u'(t), with the one-sided spectrum of P.wind.spectrum
%   von Karman  S = 4 s^2 (Lu/U) (1 + 70.8 n^2)^(-5/6),  n = f Lu/U
%   Kaimal      S = 4 s^2 (Lu/U) (1 + 6 n)^(-5/3)
% and, when P.wind.vertical is set, a correlated vertical component w'(t)
% (von Karman vertical form, sigma_w = sw_su*sigma_u, correlation r_uw).
% Each record is rescaled to zero mean and the exact target variance.
% Random phases come from a portable L'Ecuyer generator and are shared
% across mean speeds. The velocity record depends on U and the spectrum;
% it is sampled at P.wind.dt and interpolated independently of the solver step.

Wp = P.wind;
if U <= 0, error('V35:wind','U must be positive.'); end
Lu = Wp.LuL*P.L; su = Wp.Iu*U;
dt = Wp.dt; Nrec = round(Wp.T/dt) + 1;
Nfft = 2^nextpow2(Nrec);
df = 1/(Nfft*dt);
f = (1:Nfft/2-1).'*df;
nu = f*Lu/U;
switch lower(Wp.spectrum)
    case 'vonkarman'
        Su = 4*su^2*(Lu/U)*(1 + 70.8*nu.^2).^(-5/6);
    case 'kaimal'
        Su = 4*su^2*(Lu/U)./(1 + 6*nu).^(5/3);
    otherwise
        error('V35:wind','Unknown spectrum %s.', Wp.spectrum);
end
Su = Su*(su^2/(sum(Su)*df));          % discrete variance normalisation
phi = 2*pi*v35_uniform(seed, 1, numel(f));
up = synth(sqrt(2*Su*df), phi, Nfft);
up = up(1:Nrec); up = up - mean(up); up = up*(su/std(up,1));

wp = zeros(Nrec,1);
if Wp.vertical
    sw = Wp.sw_su*su; Lw = Wp.LwL*P.L; nw = f*Lw/U;
    Sw = 4*sw^2*(Lw/U)*(1 + 755.2*nw.^2)./(1 + 283.2*nw.^2).^(11/6);
    Sw = Sw*(sw^2/(sum(Sw)*df));
    gam = Wp.r_uw*su*sw/(sum(sqrt(Su.*Sw))*df);
    if abs(gam) > 1, error('V35:wind','Coherence %.3f exceeds 1.', gam); end
    psi = 2*pi*v35_uniform(seed, 2, numel(f));
    B = sqrt(2*Sw*df);
    wp = gam*synth(B, phi, Nfft) + sqrt(1-gam^2)*synth(B, psi, Nfft);
    wp = wp(1:Nrec); wp = wp - mean(wp); wp = wp*(sw/std(wp,1));
end

W = struct();
W.U = U; W.seed = seed; W.dt = dt; W.N = Nrec; W.T = (Nrec-1)*dt;
W.u = U + up; W.w = wp;
W.Ueq = v35_ueq(U, su);
W.T_ramp = Wp.T_ramp;
W.sigma_u = su; W.Iu = Wp.Iu; W.Lu = Lu; W.spectrum = Wp.spectrum;
W.vertical = Wp.vertical;
W.r_uw_realised = 0;
if Wp.vertical
    cu = W.u - mean(W.u); W.r_uw_realised = (cu.'*wp)/(norm(cu)*norm(wp));
end
end

function x = synth(A, ph, Nfft)
X = zeros(Nfft,1);
k = (2:numel(A)+1).';
X(k) = 0.5*Nfft*A.*exp(1i*ph);
X(Nfft-k+2) = conj(X(k));
x = real(ifft(X));
end
