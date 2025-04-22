%% Parameters
r      = 0.06;
sigma  = 0.3;
T      = 1;
K      = 10;
S_star = 15;

NS = 400;
Nt = 1000;

%% Solver
[S_grid, t_grid, V] = CN_PSOR(r, sigma, T, K, S_star, NS, Nt);

%% American put value V(S,t) for t = 0, 0.5, 1
figure;
hold on;
plot(S_grid, V(:,1), 'r', 'DisplayName', 't=0');
plot(S_grid, V(:,round(Nt/2)), 'g', 'DisplayName', 't=0.5');
plot(S_grid, V(:,end), 'b', 'DisplayName', 't=1');
xlabel('S'); ylabel('V(S,t)');
title('American put option');
legend show;
box on;

%% Continuation region V(S,t) > (K−S)^+
PayoffMatrix = max(K - S_grid, 0);
PayoffMat = repmat(PayoffMatrix, 1, Nt+1);
ContinuationValues = V;
ContinuationValues(V <= PayoffMat) = NaN;

[S_mat, T_mat] = meshgrid(S_grid, t_grid);

figure;
pcolor(S_mat, T_mat, ContinuationValues');
shading interp;
colormap(turbo);
colorbar('eastoutside');
xlim([7 15]);
ylim([0 1]);
xlabel('S', 'FontSize', 12);
ylabel('t', 'FontSize', 12);
title('Continuation region', 'FontWeight', 'bold');
set(gca, 'FontSize', 11);
