% Finite Element Method: 1D - Diffusion Equation 
MyU = diffusion(32,30,1,30);
close all
function U = diffusion(nx,nt,xf,tf)
clc
if nargin < 1, nx = 10; end
if nargin < 2, nt = nx; end
if nargin < 3, xf = 1; end
if nargin < 4, tf = 30; end

% spatial parameters
x0 = 0;
x = linspace(x0,xf,nx);
h = xf/(nx-1);
t0 = 0;
% time parameters
t = linspace(t0,tf,nt);
dt = tf/(nt-1);

u0 = sin(pi*x)+sin(2*pi*x);
u0(0) = 0;u0(end) = 0;

a = pi^-2; % thermal diffusivity constant


K = full(gallery('tridiag',nx,-1,2,-1))/h;
M = h*full(gallery('tridiag',nx,-1,2,-1))/6;

for i = 2:nx
    for j = 1:nt % recall nt = nx
%         U = tridiag_solve(MyL,MyD,MyU,U);
    end
end


% clf
% phi_array = zeros(1,nx+1);
% for i = 2:nx
%     i
%     phi = @(x) x+i;
%     phi = @(x) (x(i)-x)/h.*(x(i-1) < x < x(i)) + (x-x(i))/h.*(x(i) < x < x(i+1));
%     %keyboard
%     phi(x(i))
%     phi_array(i-1) = phi(x(i));
% end
% plot(x,phi_array)

keyboard

end