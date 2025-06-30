function X = rand_exponential(theta,N,seed)
    if nargin < 3, seed = 12345; end
    U = lcg_uniform(N, seed);
    X = -theta * log(U);
end