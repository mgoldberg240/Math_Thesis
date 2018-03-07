% This code was created by Matt Goldberg on Jan 12 2018
% Goal: Model the Diffusion Equations using the Crank-Nicolson method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis" (ch
% 9.1.3)


% Check analytic solution using MAPLE or by hand
% Compare to see if you get the same as Heat.pdf
% MyU = diffusion(128,128,1,1);
n = 128;
u_init = @(x) sin(pi*x);
func_U = @(x,t)(exp((-pi^2)*t)).*sin(pi*x);
diffusionCN(n,n,1,1,0.3,u_init,func_U)

function [] = diffusionCN(nt,nx,alpha,xmax,tmax,u_init,func_U)
% close all
clc
if nargin < 1, nx = 32; end
if nargin < 2, nt = 32; end
if nargin < 3, xmax = 1; end
if nargin < 4, tmax = 1; end

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

% recursion constants
a = alpha;
%a = pi^-2; % thermal diffusivity constant
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

% theoretical U
actual_U = func_U(x_vec,t_vec);
size(actual_U)
size(U)
[X,T] = meshgrid(x_vec,t_vec);

% Calculated U
figure(1)
clf
colormap(jet)
contourf(X,T,U','LineStyle','none');
xlabel('x'),ylabel('t'),title('1D Diffusion using Crank-Nicholson Method')
colorbar

% Theoretical (actual) U
figure(2)
clf
colormap(jet)
contourf(X,T,actual_U','LineStyle','none');
xlabel('x'),ylabel('t'),title('Theoretical U')
%[C,H] = contourf(X,T,StoreU);
colorbar

% error plot
figure(3)
clf
error_U = max(abs(U-actual_U));
plot(t_vec,error_U)
xlabel('x'),ylabel('time'),title('Max U error versus time')

end
