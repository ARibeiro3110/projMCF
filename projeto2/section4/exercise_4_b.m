clear; clc;

T = 1;
mu = 0.5;
sigma = 0.3;
S0 = 1;
Nsim = 5e5;
hs = 0.005*(0.5).^(0:3);

a = @(t,x) mu*x;
b = @(t,x) sigma*x;
b_x = @(t,x) sigma;
batch = 1e4;

errS_E = zeros(size(hs));  errS_M = errS_E;
errW_E = errS_E;           errW_M = errS_E;

for k = 1:numel(hs)
    h = hs(k);
	M = round(T/h);
    sumAbsE=0; sumAbsM=0; sumSE=0; sumSM=0; sumSX=0;
    sims = 0;

    while sims < Nsim
        m   = min(batch,Nsim-sims);
        Z   = randn_boxmuller(M*m, 12345+17*k+13*sims);
        Z   = reshape(Z,M,m);
        dW  = sqrt(h).*Z;

        SE  = zeros(m,1);  SM = SE;  SX = SE;
        for j=1:m
            [~,SEj] = euler_maruyama(a,b,S0,T,M,dW(:,j));
            [~,SMj] = milstein     (a,b,b_x,S0,T,M,dW(:,j));
            W       = sum(dW(:,j));
            SX(j)   = S0*exp((mu-0.5*sigma^2)*T + sigma*W);
            SE(j)   = SEj(end);   SM(j)=SMj(end);
        end

        sumAbsE = sumAbsE+sum(abs(SE-SX));
        sumAbsM = sumAbsM+sum(abs(SM-SX));
        sumSE   = sumSE+sum(SE);  sumSM=sumSM+sum(SM); sumSX=sumSX+sum(SX);
        sims    = sims+m;
    end
    errS_E(k)=sumAbsE/Nsim;   errS_M(k)=sumAbsM/Nsim;
    meanX    = sumSX/Nsim;
    errW_E(k)=abs(sumSE/Nsim-meanX);
    errW_M(k)=abs(sumSM/Nsim-meanX);
end

ord  = @(e) -diff(log(e))./log(2);
fprintf('\n   h        strong-E     strong-M      weak-E       weak-M\n');
for k=1:numel(hs)
    fprintf('%8.5f  %10.4e %10.4e %10.4e %10.4e\n', ...
          hs(k),errS_E(k),errS_M(k),errW_E(k),errW_M(k));
end
fprintf('\nEstimated orders (using last three hs):\n');
fprintf('Euler-Maruyama     strong %.3f  |  weak %.3f\n',mean(ord(errS_E)),mean(ord(errW_E)));
fprintf('Milstein  strong %.3f  |  weak %.3f\n',mean(ord(errS_M)),mean(ord(errW_M)));
