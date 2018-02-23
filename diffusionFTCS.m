% This code was created by Matt Goldberg on Jan 12 2018
% Goal: Model the Diffusion Equations using explicit (FTCS) method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis" (ch
% 9.1)
% stable = diffusion(17,257)
unstable = diffusion(33,129)

function u = diffusion(nh,nt)
clc
fprintf('nh = %f \n',nh) % number of spatial steps h
fprintf('nt = %f \n',nt) % number of time steps
a = pi^-2; % thermal diffusivity constant


x0 = 0;
xf = 1;
x = linspace(x0,xf,nh);
h = xf/(nh-1);

t0 = 0;
tf = 1; % end time
t = linspace(t0,tf,nt);
dt = tf/(nt-1);

% Stability
if dt <= h^2/(2*a)
    cond = 'stable';
    disp('h and dt satisfy stability conditions')
else
    cond = 'unstable';
    disp('h and dt DO NOT satisfy stability conditions')
end

u0 = sin(pi*x) + sin(2*pi*x);

g0 = 0;
g1 = g0; % homogenous boundary conditions

% initialize u(x,t) array
u = zeros(nh,nt);
% ICs
u(1,1:end) = g0; % u(0,t) = g0(t)
u(nh,1:end) = g1; % u(1,t) = g1(t)
u(:,1) = u0;

for n = 1:nt-1
    for i = 2:nh-1
        u(i,n+1) = u(i,n) + (a*dt/h^2) * (u(i-1,n)-2*u(i,n)+u(i+1,n));
    end
end
% u
% uf = u(:,end);
% figure
% plot(x,uf);
% keyboard

figure('name','diffusion u(x,t)','rend','painters','pos',[0 0 900 900])
hold on
clf

T = 1; % compare actual
xx = 0:0.01:1;
actual_soln = exp(-T)*sin(pi*xx)+exp(-4*T)*sin(2*pi*xx);
plot(xx,actual_soln,'k','linewidth',1) % actual solution
title(['Diffusion u(x,t) for h = ' num2str(h) ' and dt = ' num2str(dt) ' (' cond ')'])
xlabel('x'),ylabel('u')
p_c = jet(nt); % plot color

% animation
h=gca;
for i = 1:nt % we start seeing instability for h = 1/32, dt = 1/128 around t = 45 
    hold on
    h.YLimMode='manual';
    y = u(:,i);
    % axis([0 1 min(y)-5 max(y)+5])
    axis([0 1 -2 2])
    p = plot(x,y,'color',p_c(i,:));
    pause(0.001);
    set(p,'Ydata',y);
    refreshdata
    drawnow
end
% plot(x,u(:,end),'ro','linewidth',1)
end


