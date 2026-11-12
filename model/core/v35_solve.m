function x = v35_solve(K, r, damp, cV)
%V35_SOLVE Solve (K + cV*U*diag(cd)*U') x = r by the Woodbury identity.
% K is sparse and banded; the modal damping is a rank-n_damped correction.
if isempty(damp) || cV == 0
    x = K\r;
    return;
end
U = damp.U;
X = K\[r, U];
T = diag(1./(cV*damp.cd)) + U.'*X(:,2:end);
x = X(:,1) - X(:,2:end)*(T\(U.'*X(:,1)));
end
