function X = rand_uniform(a,b,N,seed)
    if nargin < 4, seed = 12345; end
    U = lcg_uniform(N, seed);
    X = a + (b-a) .* U;
end
