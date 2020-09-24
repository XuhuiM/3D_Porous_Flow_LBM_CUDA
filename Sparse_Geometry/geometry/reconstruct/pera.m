clear;

load u.dat;

load v.dat;

load rho.dat;

% load flag.dat;

% streamslice(u ,v, 10);
% 
% daspect([1 1 1]);

L = 1.0;

N = 200;

% Nc = 60;
% % 
% R = 0.1;

dx = L / N;

Dp = 0.4;

tau_f = 0.8;

niu = (tau_f - 0.5) * dx / 3.0;

% Red = 0.1;
% 
% u0 = Red * niu / Dp;
% 
% phi = 0.25 * pi * (Dp / L)^2
% 
% epsilon = 1.0 - phi;
% 
% ka = Dp * Dp * (-log(phi) - 1.476 + 2.0 * phi - 1.774 * phi * phi + 4.706 * phi * phi * phi) / 32 / phi;
% 
% k = Dp * Dp * (-log(phi) - 1.476 + 2.0 * phi - 1.774 * phi * phi + 4.706 * phi * phi * phi);
% 
% G = niu * u0 / k;

G = 1.0e-6;

% [col len] = size(rho);
% 
% for i = 1 : len
%     
%     sum_rho = 0.0;
%     
%     sum_n = 0;
%     
%     for j = 1 : col
%         
%         if (flag(j,i) == 1)
%           
%             sum_rho = sum_rho + rho(j, i);
%             
%             sum_n = sum_n + 1;
%                        
%         end
%         
%         arho(i) = sum_rho / sum_n;
%                
%     end
%     
% end

% k = 1 : len;
% 
% x = (k - 1) * dx;
% 
% a = polyfit(x, arho / 3.0, 1);
% 
% % G = - a(1);
% % 
% % yy = polyval(a, x);
% 
% rho_in  = arho(1);
% 
% rho_out = arho(N+1);
% % 
% G = (rho_in - rho_out) / L / 3.0;
% 
% % plot(x, arho / 3.0, 'b-', 'linewidth', 2.0);
% % 
% % hold on;
% % 
% % plot(x, yy, 'r--', 'linewidth', 2.0);

[col len] = size(flag);

sum_n = (N+1) * (N+1); 

% for j = 1 : col
%     for i = 1 : len
%         
%         if(flag(j,i) == 0)
%             
% %             sum = sum + 1;
%             
%             sum_u = sum_u + u(j,i);
%         end
%     end
% end

sum_u = sum(sum(u));
au = sum_u / sum_n;

k = niu * au / G;

fprintf('permeability is %e\n', k);

% phi = 0.25 * pi * (Dp)^2
% 
% epsilon = 1.0 - phi;

% ka = phi * phi * phi * Dp * Dp / 150 / (1.0 - phi)^2;

% ka = Dp * Dp * (-log(phi) - 1.476 + 2.0 * phi - 1.774 * phi * phi + 4.706 * phi * phi * phi) / 32 / phi;

% ka = epsilon * epsilon * epsilon * (epsilon - 0.2146) * Dp * Dp / 31 / phi^1.3;

% fprintf('Analytical permeability is %e\n', ka);

% err = abs(k - ka)/ka*100