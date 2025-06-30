function U = lcg_uniform(N, seed)
    M = 2^31 - 1;
    a = 16807;
    b = 0;

    m = zeros(N,1);
    m(1) = mod(seed, M); % garante 0 < m1 < M

    for k = 2:N
        m(k) = mod(a * m(k-1) + b, M);
    end

    U = m / M;
end
