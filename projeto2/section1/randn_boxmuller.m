function Z = randn_boxmuller(N, seed)
    if nargin < 2, seed = 12345; end
    U = lcg_uniform(N, seed);
    U1 = U(1:floor(N/2));
    U2 = U(floor(N/2)+1:end);

    R     = sqrt(-2 .* log(U1));
    Theta = 2*pi .* U2;

    Z = [R .* cos(Theta); R .* sin(Theta)];

end
