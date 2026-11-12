function [lam, Phi] = v35_eig_low(K, M, k)
%V35_EIG_LOW Smallest k eigenvalues of K phi = lam M phi (real parts, ascending).
nf = size(K,1); k = min(k, nf);
if nf <= 120
    [V, D] = eig(full(K), full(M));
    d = real(diag(D));
    [d, ix] = sort(d); lam = d(1:k); Phi = real(V(:, ix(1:k)));
else
    sig = 0;                          % shift-invert about zero
    opts.tol = 1e-12; opts.maxit = 500; opts.disp = 0;
    kk = min(nf-2, k+4);
    [V, D] = eigs(K, M, kk, sig, opts);
    d = real(diag(D));
    [d, ix] = sort(d); lam = d(1:k); Phi = real(V(:, ix(1:k)));
end
end
