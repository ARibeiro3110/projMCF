function H = halton2d(N)
    base = [2,3];
    H = zeros(N,2);
    for i = 1:N
        n = i;
        for d = 1:2
            f = 1/base(d);
            x = 0;
            while n > 0
                x = x + mod(n, base(d)) * f;
                n = floor(n/base(d));
                f = f/base(d);
            end
            H(i,d) = x;
            n = i;
        end
    end
end
