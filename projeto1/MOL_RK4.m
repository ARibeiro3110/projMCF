function out = MOL_RK4(r, sigma, T, S_star, NS, Nt, u_0, u_a, u_b)
    %% Spatial grid
    hS  = S_star / NS;
    h   = T / Nt;
    S   = 0 : hS : S_star;
    i   = (1:NS-1)';


    %% Finite-difference coefficients α_i, β_i, γ_i
    alpha_i = 0.5*sigma^2 .* i.^2  - 0.5*r .* i;
    beta_i  =    -sigma^2 .* i.^2  -     r;
    gamma_i = 0.5*sigma^2 .* i.^2  + 0.5*r .* i;

    A = spdiags([alpha_i beta_i gamma_i], [-1 0 1], NS-1, NS-1); % A_ML


    %% Initial condition W(0) (pay-off)
    W = u_0(S(2:NS)).';


    %% Pre-allocate storage for all slices
    U = zeros(NS+1, Nt+1);
    U(:,1) = u_0(S).'; % (t=0 forward)
    b = zeros(NS-1,1); % b_ML(t)


    %% RK-4 loop (forward time)
    for n = 0:Nt-1
        t = n * h;

        % Boundary values V(S=0,t) and V(S=S*,t)
        V_left  = u_a(t); % U(0,t)
        V_right = u_b(t); % U(S*,t)
        b(1) = alpha_i(1) * V_left;

        % Runge-Kutta stages
        f1 = h * (A*W + b);

        t_half   = t + 0.5*h;
        V_left_h = u_a(t_half);
        b(1)     = alpha_i(1) * V_left_h;
        f2 = h * (A*(W + 0.5*f1) + b);
        f3 = h * (A*(W + 0.5*f2) + b);

        t_next   = t + h;
        V_left_n = u_a(t_next);
        b(1)     = alpha_i(1) * V_left_n;
        f4 = h * (A*(W + f3) + b);

        % RK-4 update
        W = W + (f1 + 2*f2 + 2*f3 + f4)/6;

        % Assemble full vector V(S_i, t_next) and store
        U(:,n+2) = [V_left_n ; W ; u_b(t_next)];
    end


    %% Output
    out.S = S;
    out.t = 0:h:T;
    out.U = U;
end
