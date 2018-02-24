% This code was created by Matt Goldberg on Jan 12 2018
% Goal: Model the Diffusion Equations using the Crank-Nicolson method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis" (ch
% 9.1.3)


% Check analytic solution using MAPLE or by hand
% Compare to see if you get the same as Heat.pdf
MyU = diffusion(30,30,1,1);

function U = diffusion(nx,nt,xf,tf)
close all
clc
if nargin < 1, nx = 10; end
if nargin < 2, nt = 30; end
if nargin < 3, xf = 1; end
if nargin < 4, tf = 30; end

fprintf('nx = %f \n',nx)
fprintf('nt = %f \n',nt)

x0 = 0;
x_vec = linspace(x0,xf,nx);
h = xf/nx;

t0 = 0;
t_vec = linspace(t0,tf,nt);
dt = tf/nt;
u0 = sin(pi*x_vec)+sin(2*pi*x_vec);


a = pi^-2; % thermal diffusivity constant
r = a*dt/(2*(h^2));

MyLower = -r*ones(1,nx); MyLower(1) = 0;
MyDiag = (2*r+1)*ones(1,nx); MyDiag(1) = 1; MyDiag(end) = 1;
MyUpper = -r*ones(1,nx); MyUpper(end) = 0;

U(:,1) = u0; 

for t = 2:nt
    W1 = -[0 MyLower(2:end-1).*U(1:end-2,t-1)' 0];
    W2 = [0 MyDiag(2:end-1).*U(2:end-1,t-1)' 0];
    W3 = -[0 MyUpper(2:end-1).*U(3:end,t-1)' 0];
    W = W1+W2+W3;

    U(:,t) = tridiag_solve(MyLower,MyDiag,MyUpper,W)
    % U(:,t) = tridiag_solve(MyLower,MyDiag,MyUpper,U(:,t-1));
end

% theoretical U
actual_U = (exp(t_vec)).*sin(pi*x_vec)' + (exp(-4*t_vec)).*sin(2*pi*x_vec)';
size(actual_U)
size(U)
[X,T] = meshgrid(x_vec,t_vec);

figure(1)
clf
colormap(jet)
contourf(X,T,U','LineStyle','none');
xlabel('x'),ylabel('t'),title('1D Diffusion using Crank-Nicholson Method')
%[C,H] = contourf(X,T,StoreU);
colorbar

figure(2)
clf
colormap(jet)
contourf(X,T,actual_U','LineStyle','none');
xlabel('x'),ylabel('t'),title('1D Diffusion using Crank-Nicholson Method')
%[C,H] = contourf(X,T,StoreU);
colorbar

% error plot
figure(3)
clf
error_U = max(abs(U-actual_U));
plot(t_vec,error_U)
xlabel('x'),ylabel('time'),title('Max U error versus time')

end
