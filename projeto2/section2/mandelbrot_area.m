clear; clc;

%% (i) Load Mandelbrot
D   = load('mandelbrot.mat');
fn  = fieldnames(D);
M   = D.(fn{1});
nGrid = size(M,1);

%% (ii) visualize fractal
figure('Color','w');
imagesc([0 1],[0 1],flipud(M'));
set(gca,'YDir','normal');
axis image;
colormap(flipud(gray)); clim([0 1]);
grid on;

%% (iii) parameters
N      = 1e5;
seedMC = 54321;

%% (iv) Monte Carlo - uniform points via LCG
U = lcg_uniform(2*N, seedMC);
P = reshape(U, [N 2]);
chiMC   = chi_mandelbrot(P(:,1), P(:,2), M);
areaMC  = mean(chiMC);
stderrMC = sqrt(var(chiMC)/N);

%% (v) Quasi-Monte Carlo - Halton sequence
H = halton2d(N);
chiQMC  = chi_mandelbrot(H(:,1), H(:,2), M);
areaQMC = mean(chiQMC);

%% (vi) Results in the prompt
fprintf('--- Estimate of the Mandelbrot area ---\n');
fprintf('Monte Carlo Estimate  : %.6f   (± %.6f)\n', areaMC, stderrMC);
fprintf('Quasi-Monte Carlo Estimate (Halton): %.6f\n', areaQMC);
