function [S_grid, t_grid, V] = CN_PSOR(r, sigma, T, K, S_star, NS, Nt)
    %% Grids
    dS     = S_star / NS;
    dt     = T / Nt;
    S_grid = linspace(0, S_star, NS+1)';
    t_grid = linspace(0, T, Nt+1);

    %% Grid indices
    I = 2:NS;
    S = S_grid(I);
    M = NS - 1;

    %% Coefficients for tridiagonal matrices
    i = (1:M)';  % i indexes interior points
    alpha = 0.25 * dt * (sigma^2 * i.^2 - r * i);
    beta  = -0.5 * dt * (sigma^2 * i.^2 + r);
    gamma = 0.25 * dt * (sigma^2 * i.^2 + r * i);

    main_diag = 1 - beta;
    lower_diag = -alpha(2:end);
    upper_diag = -gamma(1:end-1);

    %% Matrix B (I - dt/2*A)
    B = spdiags([[lower_diag; 0], main_diag, [0; upper_diag]], -1:1, M, M);

    %% Matrix A+ (I + dt/2*A)
    main_plus = 1 + beta;
    low_plus = alpha(2:end);
    up_plus  = gamma(1:end-1);
    A_plus = spdiags([[low_plus; 0], main_plus, [0; up_plus]], -1:1, M, M);

    %% Payoff and boundary conditions
    V = zeros(NS+1, Nt+1);
    V(:, end) = max(K - S_grid, 0);
    V(1, :) = K; V(end, :) = 0;
    payoff = max(K - S, 0);

    %% Relaxation
    omega = 1.3;
    tol = 1e-7;
    max_iter = 500;

    %% Time-stepping loop
    for n = Nt:-1:1
        rhs = A_plus * V(I, n+1);

        % Boundary conditions
        rhs(1)   = rhs(1)   + alpha(1)   * (V(1, n+1) + V(1, n));
        rhs(end) = rhs(end) + gamma(end) * (V(end, n+1) + V(end, n));

        % Initial guess = last solution
        x = V(I, n+1);

        % PSOR iterations
        for it = 1:max_iter
            x_old = x;
            for i = 1:M
                if i == 1
                    residual = rhs(i) - B(i,i)*x(i) - B(i,i+1)*x(i+1);
                elseif i == M
                    residual = rhs(i) - B(i,i-1)*x(i-1) - B(i,i)*x(i);
                else
                    residual = rhs(i) - B(i,i-1)*x(i-1) - B(i,i)*x(i) - B(i,i+1)*x(i+1);
                end
                x(i) = max(payoff(i), x(i) + omega * residual / B(i,i));
            end
            if norm(x - x_old, inf) < tol
                break;
            end
        end

        V(I, n) = x;
    end
end
