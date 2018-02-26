% This code was created by Matt Goldberg on Jan 12 2018
% Goal: Model the Diffusion Equations using the Crank-Nicolson method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis" (ch
% 9.1.3)


% Check analytic solution using MAPLE or by hand
% Compare to see if you get the same as Heat.pdf
MyU = diffusion(128,128,1,1);

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
h = xf/(nx-1);

t0 = 0;
t_vec = linspace(t0,tf,nt);
dt = tf/(nt-1);
u0 = sin(pi*x_vec)+sin(2*pi*x_vec);


a = pi^-2; % thermal diffusivity constant
r = a*dt/(2*(h^2));

K = full(gallery('tridiag',nx,-1,2,-1))/h;

LHS_Lower = -r*ones(1,nx); LHS_Lower(1) = 0;
LHS_Diag = (2*r+1)*ones(1,nx); %LHS_Diag(1) = 1; LHS_Diag(end) = 1;
LHS_Upper = -r*ones(1,nx); LHS_Upper(end) = 0;

U(:,1) = u0; 

for t = 2:nt
%     W1 = -[0 MyLower(2:end-1).*U(1:end-2,t-1)' 0];
%     W2 = [0 MyDiag(2:end-1).*U(2:end-1,t-1)' 0];
%     W3 = -[0 MyUpper(2:end-1).*U(3:end,t-1)' 0];
%     W = W1+W2+W3;
%     
    % LHS = ((eye(nx)-r*K)^-1)*(eye(nx)+r*K);

    Un = U(:,t-1); % vector in x
    RHS_Lower = [0 -LHS_Lower(2:end-1).*Un(1:end-2)' 0];
    % RHS_Lower = Un(1) + Un( LHS_Lower;RHS_Lower(1) = 0;RHS_Lower(end) = gn;
    k = (1-2*r)/(2*r+1);
    RHS_Diag = [0 k*LHS_Diag(2:end-1).*Un(2:end-1)' 0];
    RHS_Upper = [0 -LHS_Lower(2:end-1).*Un(3:end)' 0];
    %RHS_Diag = U(:,t-1) - (1-2*r)*ones(1,nx);RHS_Diag(1) = 0;RHS_Diag(end) = 0;
    %RHS_Upper = U(:,t-1) - LHS_Lower;RHS_Lower(1) = 0;RHS_Lower(end) = 0;
    (RHS_Lower+RHS_Diag+RHS_Upper);
    RHS = Un-(RHS_Lower+RHS_Diag+RHS_Upper);
    % U(:,t) = linsolve(LHS,U(:,t-1));
    % U(:,t) = tridiag_solve(MyLower,MyDiag,MyUpper,W)
    
    %Solve for U_n+1
    U(:,t) = tridiag_solve(LHS_Lower,LHS_Diag,LHS_Upper,RHS);
    % U(:,t) = tridiag_solve(LHS_Diag,LHS_Lower,LHS_Upper,RHS);
end

% theoretical U
actual_U = (exp(-t_vec)).*sin(pi*x_vec)' + (exp(-4*t_vec)).*sin(2*pi*x_vec)';
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
