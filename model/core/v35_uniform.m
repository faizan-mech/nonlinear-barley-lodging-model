function r = v35_uniform(seed, stream, N)
%V35_UNIFORM N uniform numbers in (0,1) from seed and stream.
% Combined multiplicative congruential generator (L'Ecuyer 1988), exact in
% double precision, so MATLAB and Octave give identical sequences and the
% global random state is never touched.
m1 = 2147483563; m2 = 2147483399;
s1 = mod(seed*7919 + stream*104729 + 12345, m1-1) + 1;
s2 = mod(seed*4099 + stream*65537 + 67890, m2-1) + 1;
for k = 1:64
    s1 = mod(40014*s1, m1); s2 = mod(40692*s2, m2);
end
r = zeros(N,1);
for k = 1:N
    s1 = mod(40014*s1, m1); s2 = mod(40692*s2, m2);
    z = s1 - s2; if z < 1, z = z + m1 - 1; end
    r(k) = z/m1;
end
end
