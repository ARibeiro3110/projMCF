clear; clc; close all;

T = 1;
h = 1e-3;
N = T/h;
mu = 0.5;
sigma = 0.3;
S0 = 1;

a = @(t,x) mu*x;
b = @(t,x) sigma*x;
b_x = @(t,x) sigma;

dW = sqrt(h) * randn_boxmuller(N,12345);

[ t, SE ] = euler_maruyama(a,b,S0,T,N,dW);
[ ~, SM ] = milstein(a,b,b_x,S0,T,N,dW);

W      = cumsum([0; dW]);                   % W_0=0
SExact = S0*exp((mu-0.5*sigma^2)*t + sigma*W);

plot(t,SExact,'k', t,SE,'b--', t,SM,'r:','LineWidth',1.2);
xlabel('t'); ylabel('S(t)');
legend('Exact','Euler','Milstein');
title('GBM: exact vs numerical');
grid on;
