% This code was created by Matt Goldberg on Jan 12 2018
% Goal: Model the Diffusion Equations using explicit (FTCS) method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis" (ch
% 9.1)
function u = diffusionFTCS(h,dt)
clc
fprintf('h = %f \n',h)
fprintf('dt = %f \n',dt)
a = pi^-2; % thermal diffusivity constant

% Stability
if dt <= h^2/(2*a)
    cond = 'stable';
    disp('h and dt satisfy stability conditions')
else
    cond = 'unstable';
    disp('h and dt DO NOT satisfy stability conditions')
end

% h = 1/32; % spatial step size
x0 = 0;
xf = 1;
x = x0:h:xf;
x_len = length(x); % used for iteration
% x = linspace(x0,xf,x_len-1); % spatial vector

% dt = 1/128;
t0 = 0;
tf = 1; % end time
t = t0:dt:tf;
t_len = tf/dt;
% t = linspace(t0,tf,t_len); % time vector

u0 = sin(pi*x) + sin(2*pi*x);

g0 = 0;
g1 = g0; % homogenous boundary conditions

% initialize u(x,t) array
u = zeros(length(x),t_len);

% ICs
u(1,1:end) = g0; % u(0,t) = g0(t)
u(x_len,1:end) = g1; % u(1,t) = g1(t)
u(:,1) = u0;


x_ind = x0 + 2; % index = 2 --- 1st row is filled
t_ind = t0 + 1; % index = 1 --- 1st column

for n = t_ind:t_len-1
    %i
    for i = x_ind:x_len-1
        %n
        %clc
        u(i,n+1) = u(i,n) + (a*dt/h^2) * (u(i-1,n)-2*u(i,n)+u(i+1,n));
        %pause(0.5)
    end
end



figure('name','diffusion u(x,t)','rend','painters','pos',[0 0 900 900])
hold on
clf

T = 1; % compare actual
xx = 0:0.01:1;
actual_soln = exp(-T)*sin(pi*xx)+exp(-4*T)*sin(2*pi*xx);
plot(xx,actual_soln,'k','linewidth',1) % actual solution
title(['Diffusion u(x,t) for h = ' num2str(h) ' and dt = ' num2str(dt) ' (' cond ')'])
xlabel('x'),ylabel('u')
p_c = jet(t_len); % plot color

% animation
h=gca;
for i = 1:t_len % we start seeing instability for h = 1/32, dt = 1/128 around t = 45 
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


