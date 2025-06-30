function [t_grid,X] = euler_maruyama(a,b,X0,T,N,dW)
    if nargin < 6
        dW = sqrt(T/N)*randn_boxmuller(N);
    end

    h      = T/N;
    t_grid = linspace(0,T,N+1).';
    X      = zeros(N+1,1);   X(1)=X0;
    
    for n = 1:N
        t        = t_grid(n);
        X(n+1)   = X(n) + a(t,X(n))*h + b(t,X(n))*dW(n);
    end
end
