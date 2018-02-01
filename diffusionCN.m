% This code was created by Matt Goldberg on Jan 12 2018
% Goal: Model the Diffusion Equations using the Crank-Nicolson method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis" (ch
% 9.1.3)
function u = diffusionCN(nx,nt)
clc
if nargin < 1, nx = 50; end
if nargin < 2, nt = 50; end   

fprintf('nx = %f \n',nx)
fprintf('nt = %f \n',nt)
a = 0.1; % thermal diffusivity constant

x0 = 0;
xf = 1;
x = linspace(x0,xf,nx);
h = xf/(nx-1);

t0 = 0;
tf = 30; % end time
t = linspace(t0,tf,nt);
dt = tf/(nt-1);

u0 = 4.*x.*(1-x);

% homogenous boundary conditions
g0 = 0; u(1) = g0;
g1 = g0; u(end) = g1;
a = 1; % thermal diffusivity constant
r = a*dt/(nx^2);

MyL = -r*ones(1,nx); MyL(1) = 0;
MyD = (2*r+1)*ones(1,nx);
MyU = -r*ones(1,nx); MyU(end) = 0;
u = tridiag_solve(MyL,MyD,MyU,u0); % calculate u1 from u0


U = [];
U(:,1) = u;

for n = 2:nt
    U(:,n) = tridiag_solve(MyL,MyD,MyU,U(:,n-1)); 
end

U = U';
[X T] = meshgrid(x,t);
figure
size(X),size(T),size(U)
contour(X,T,U)
colorbar
