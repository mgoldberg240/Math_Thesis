% This code was created by Matt Goldberg on Jan 12 2018
% Goal: Model the Diffusion Equations using the Crank-Nicolson method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis" (ch
% 9.1.3)


% Check analytic solution using MAPLE or by hand
% Compare to see if you get the same as Heat.pdf
MyU = diffusion(10,30,1,30);

function U = diffusion(nx,nt,xf,tf)
clc
if nargin < 1, nx = 10; end
if nargin < 2, nt = 30; end
if nargin < 3, xf = 1; end
if nargin < 4, tf = 30; end

fprintf('nx = %f \n',nx)
fprintf('nt = %f \n',nt)

x0 = 0;
% xf = 1;
x = linspace(x0,xf,nx);
h = xf/(nx-1);

t0 = 0;
% tf = 15; % end time
t = linspace(t0,tf,nt);
dt = tf/(nt-1);
u0 = 4*x.*(1-x);%2-1.5.*x+sin(pi*x); % sin(pi*x);


a = 1;%pi^-2; % thermal diffusivity constant
r = a*dt/(nx^2);

MyL = -r*ones(1,nx); MyL(1) = 0;
MyD = (2*r+1)*ones(1,nx); MyD(1) = 1; MyD(end) = 1;
MyU = -r*ones(1,nx); MyU(end) = 0;
u = tridiag_solve(MyL,MyD,MyU,u0); % calculate u1 from u0


U(1,:) = u; 
StoreU = U;
for n = 2:nt
%     W1 = -[0 MyL(2:end-1).*U(1:end-2) 0];
%     W2 = [0 2*(1-MyL(2:end-1)).*U(2:end-1) 0];
%     W3 = - [0 MyU(2:end-1).*U(3:end) 0];
%     W = W1+W2+W3;
%     MyD(1) = g0; MyD(end) = g1;
    U = tridiag_solve(MyL,MyD,MyU,U);
    
    StoreU(n,:) = U;
end

[X,T] = meshgrid(x,t);
figure
clf
colormap(jet)
contourf(X,T,StoreU,'LineStyle','none');
xlabel('x'),ylabel('t'),title('1D Diffusion using Crank-Nicholson Method')
%[C,H] = contourf(X,T,StoreU);
colorbar
end
