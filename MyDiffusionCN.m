% This code was created by Matt Goldberg on Jan 12 2018
% Last updated Mar 6 2018
% Goal: Model the Diffusion Equations using the Finite Element method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis"
% Similar to Crank Nicholson method (MyDiffusionCN)

% n = 128;
% u_init = @(x) sin(pi*x);
% % func_U = @(x,t)(exp((-pi^2)*t)).*sin(pi*x);
% % diffusionCN(n,n,1,1,0.3,u_init,func_U)

function [U,actual_U,E,X,T] = MyDiffusionCN(nt,nx,a,xmax,tmax,u_init,func_U)
% close all
clc

% --- x grid
x0 = 0;
x_vec = linspace(x0,xmax,nx)';
h = xmax/(nx-1);
% --- t grid
t0 = 0;
t_vec = linspace(t0,tmax,nt);
dt = tmax/(nt-1);

% recursion constants
% a = alpha = thermal diffusivity constant 

% --- coefficient in tridiag system
r = a*dt/(2*(h^2));

LHS_Lower = -r*ones(nx,1); LHS_Lower(end) = 0;
LHS_Diag = (2*r+1)*ones(nx,1); LHS_Diag(1) = 1; LHS_Diag(end) = 1;
LHS_Upper = -r*ones(nx,1); LHS_Upper(1) = 0;


% --- u0 
U(:,1) = u_init(x_vec); 
k = (1-2*r)/(2*r+1);

for t = 2:nt
    % current U_n
    Un = U(:,t-1);
    
    % construct Right Hand Side vector corresponding to t_n
    RHS_Lower = [0; -LHS_Lower(2:end-1).*Un(1:end-2); 0];
    RHS_Diag = [0; k*LHS_Diag(2:end-1).*Un(2:end-1); 0];
    RHS_Upper = [0; -LHS_Lower(2:end-1).*Un(3:end); 0];
    RHS = (RHS_Lower+RHS_Diag+RHS_Upper);
    
    % set boundary conditions
    RHS(1) = 0; RHS(end) = 0;

    %Solve for U_n+1
    U(:,t) = tridiag_solve(LHS_Lower,LHS_Diag,LHS_Upper,RHS);
end
U = U';

% theoretical U
actual_U = func_U(x_vec,t_vec);
[X,T] = meshgrid(x_vec,t_vec);
E = max(abs(U'-actual_U)); % error

% % Calculated U
% figure(1)
% clf
% colormap(jet)
% AC = contourf(X,T,U,'LineStyle','none');
% xlabel('x'),ylabel('t'),title('1D Diffusion using Crank-Nicholson Method')
% colorbar
% 
% % Theoretical (actual) U
% figure(2)
% clf
% colormap(jet)
% TC = contourf(X,T,actual_U','LineStyle','none');
% xlabel('x'),ylabel('t'),title('Theoretical U')
% %[C,H] = contourf(X,T,StoreU);
% colorbar
% 
% % error plot
% figure(3)
% clf
% E = max(abs(U'-actual_U)); % error
% plot(t_vec,E)
% xlabel('x'),ylabel('time'),title('Max U error versus time')

end
