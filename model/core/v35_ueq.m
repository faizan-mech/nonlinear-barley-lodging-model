function Ueq = v35_ueq(U, sigma)
%V35_UEQ Speed with Ueq^2 = <u|u|> for Gaussian u ~ N(U, sigma^2).
% The mean-wind equilibrium uses Ueq on the ear and alpha(s)*Ueq on the culm,
% so the equilibrium carries the mean of the quadratic load.
if sigma <= 0, Ueq = U; return; end
a = U/sigma;
Phi = 0.5*(1 + erf(a/sqrt(2)));
phi = exp(-0.5*a^2)/sqrt(2*pi);
Ueq = sqrt(max((U^2 + sigma^2)*(2*Phi - 1) + 2*U*sigma*phi, 0));
end
