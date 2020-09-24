clear;

load u.dat;
% streamslice(u ,v, 10);
% 
% daspect([1 1 1]);

L = 1.0;

N = 90;

% Nc = 60;
% % 
% R = 0.1;

dx = L/N;

tau_f = 1.0;

niu = (tau_f - 0.5) * dx / 3.0;

sum_n = (N+1)^3; 

sum_u = 0.0;

sum_u = sum(sum(u));

au = sum_u / sum_n;

G = 1.0e-5;

k = niu * au / G;

fprintf('permeability is %e\n', k);