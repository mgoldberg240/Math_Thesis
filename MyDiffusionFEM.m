% This code was created by Matt Goldberg on Jan 12 2018
% Last updated Apr 2 2018
% Goal: Model the Diffusion Equations using the Crank-Nicolson method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis" (ch
% 9.1.3)


% Check analytic solution using MAPLE or by hand
% Compare to see if you get the same as Heat.pdf
% 
% n = 33;
% u_init = @(x) sin(pi*x)+sin(2*pi*x);
% func_U = @(x,t)(exp(-t)).*sin(pi*x) + (exp(-4*t)).*sin(2*pi*x);
% [U_FEM,E_FEM] = DiffusionFEM(n,n,pi^-2,1,1,u_init,func_U);

%function [U,E] = DiffusionFEM(nt,nx,a,xmax,tmax,u_init,func_U)

function [U,E] = MyDiffusionFEM(nt,nx,a,xmax,tmax,u_init,func_U)
% close all

% if nargin < 1, nx = 32; end
% if nargin < 2, nt = 32; end
% if nargin < 3, xmax = 1; end
% if nargin < 4, tmax = 1; end

% fprintf('nx = %f \n',nx)
% fprintf('nt = %f \n',nt)

% --- x grid
x0 = 0;
x_vec = linspace(x0,xmax,nx)';
h = xmax/(nx-1);
% --- t grid
t0 = 0;
t_vec = linspace(t0,tmax,nt);
dt = tmax/(nt-1);

%a = pi^-2; % thermal diffusivity constant

% --- recursion coefficients in tridiag system
r = a*dt/2;
q = h/6;
s = r/h;


LHS_Lower = (q-s)*ones(nx,1); LHS_Lower(end) = 0;
LHS_Diag = (4*q+2*s)*ones(nx,1); LHS_Diag(1) = 1; LHS_Diag(end) = 1;
LHS_Upper = (q-s)*ones(nx,1); LHS_Upper(1) = 0;


% --- u0 
U(:,1) = u_init(x_vec); 
k1 = (q+s)/(q-s);
k2 = (4*q-2*s)/(4*q+2*s);

for t = 2:nt
    % current U_n
    Un = U(:,t-1);
    
    % construct Right Hand Side vector corresponding to t_n
    RHS_Lower = [0; k1*LHS_Lower(2:end-1).*Un(1:end-2); 0];
    RHS_Diag = [0; k2*LHS_Diag(2:end-1).*Un(2:end-1); 0];
    RHS_Upper = [0; k1*LHS_Lower(2:end-1).*Un(3:end); 0];
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

% figure
% clf
% hold on
% plot(t_vec,E)

% figure
% surf(X,T,U')
% figure
% hold on
% plot(t_vec,E,'k')
end
