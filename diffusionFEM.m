% Finite Element Method: 1D - Diffusion Equation 
MyU = diffusion(32,32,1,1);

function U = diffusion(nx,nt,xf,tf)
close all
clc
if nargin < 1, nx = 10; end
if nargin < 2, nt = nx; end
if nargin < 3, xf = 1; end
if nargin < 4, tf = 30; end

% spatial parameters
x0 = 0;
x_vec = linspace(x0,xf,nx);
h = xf/nx;
t0 = 0;
% time parameters
t_vec = linspace(t0,tf,nx); % currently nx is size as nt -> square result U
dt = tf/nt-1;

u0 = sin(pi*x_vec)+sin(2*pi*x_vec);
u0(1) = 0;
u0(end) = 0;
U = zeros(nx);
U(:,1) = u0;
a = pi^-2; % thermal diffusivity constant


K = full(gallery('tridiag',nx,-1,2,-1))/h;
k = a*dt/2;
M = h*full(gallery('tridiag',nx,1,4,1))/6;



for t = 2:nt
%         % LHS = (M-k*K)^-1*(M+k*K);
        LHS = eye(nx)+dt*M^-1*K^-1;
%         U = tridiag_solve(MyL,MyD,MyU,U);

        % solve LHS*U_i+1 = U_i
        U(:,t) = linsolve(LHS,U(:,t-1));
end

[X,T] = meshgrid(x_vec,t_vec);
% size(X),size(T),size(U)

% theoretical U

actual_U = (exp(t_vec)).*sin(pi*x_vec)' + (exp(-4*t_vec)).*sin(2*pi*x_vec)';
size(actual_U)
figure(1)
colormap(jet)
contourf(X,T,actual_U,'LineStyle','none');
xlabel('x'),ylabel('time'),title('Theoretical U(x,t)')
colorbar

figure(2)
contourf(X,T,U,'LineStyle','none');
xlabel('x'),ylabel('time'),title('Calculated U(x,t)')
% [C,H] = contourf(X,T,U');
colormap(jet)
colorbar

% error plot
figure(3)
error_U = max(abs(U-actual_U));
plot(t_vec,error_U)
xlabel('x'),ylabel('time'),title('Max U error versus time')
keyboard

end