function [t_grid,X] = milstein(a,b,b_deriv,X0,T,N,dW)
    if nargin < 7
        dW = sqrt(T/N)*randn_boxmuller(N);
    end

    h      = T/N;
    t_grid = linspace(0,T,N+1).';
    X      = zeros(N+1,1);   X(1)=X0;
    
    for n = 1:N
        t      = t_grid(n);
        x      = X(n);
        dw     = dW(n);
        X(n+1) = x + a(t,x)*h + b(t,x)*dw ...
                   + 0.5*b(t,x)*b_deriv(t,x)*(dw^2-h);
    end
end
