%% Parameters
r      = 0.06;
sigma  = 0.3;
T      = 1;
K      = 10;
S_star = 15;

NS = 400;
Nt = 13000;

%% Functions for European put option
u_0 = @(S) max(K-S,0);
u_a = @(t) K*exp(-r*t);
u_b = @(t) 0*t;

%% Run
sol = MOL_RK4(r,sigma,T,S_star,NS,Nt,u_0,u_a,u_b);

%% 2D price at calendar time t=0 (last column)
V_today = sol.U(:,end);
figure;
plot(sol.S, V_today), grid on
title('European put V(S,0)')

%% 3D surface V(S,t)
[Tgrid, Sgrid] = meshgrid(T - sol.t, sol.S); % calendar time axis
figure;
mesh(Tgrid, Sgrid, sol.U)
xlabel('t'),
ylabel('S'),
zlabel('V(S,t)')
title('European put option value, V(S,t)')
view(135,30);
