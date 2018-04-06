% This code was created by Matt Goldberg on Jan 12 2018
% Goal: Model the Diffusion Equations using explicit (FTCS) method as presented
% By Epperson in "An Introduction to Numerical Methods and Analysis" (ch
% 9.1)

u_init = @(x) sin(pi*x)+sin(2*pi*x); % u0
func_U = @(x,t)(exp(-t)).*sin(pi*x) + (exp(-4*t)).*sin(2*pi*x); % actual u(x,t)


%% Plot error for changing nt
close all
figure('name','FTCS errors','rend','painters','pos',[0 0 1000 1000]);
clf; hold on
g=gca;
set(gca,'linewidth',3,'fontsize',20)
g.TickLength = [0.02 0.02];


Nplots = 4;
Stab = cell(1,Nplots);
Error = zeros(1,Nplots);
nx = 32;
nt = 256;
nx_vals = nx*ones(1,Nplots);
nt_vals = 16*[16 64 256 1024];
[U,E,cond,t] = DiffusionFTCS(nx,nt,a,1,1,u_init,func_U);
style = {'.';'+';'--';'o'};
% colors = {'r','g','b','m'};
% colors = [46 53 50; 126 145 129; 199 206 219; 123,211,137; 148 132 155];
colors = jet(Nplots);
colors(3,:) = [100 255 0]/255;
for i = 1:Nplots
    a = pi^-2;
    nx = nx_vals(i);
    nt = nt_vals(i);%floor(nx^2/(2*a));
    
    [U,E,cond,t] = DiffusionFTCS(nx,nt,a,1,1,u_init,func_U);

    Stab{i} = ['$h = 1/$' num2str(nx) ', $\Delta t = 1/$' num2str(nt) ' --- ' num2str(cond)];
    Error(i) = max(E(:));
    p(i) = plot(t,E,'color',colors(i,:),'LineStyle','-','LineWidth',5);
end

% leg_str{i} = ['nx = ' num2str(nx_vals(i)) ', nt = ' num2str(nx_vals(i))];
lgd1 = legend(Stab);
set(lgd1,'Interpreter','LaTeX');
lgd1.FontSize=22;

t1 = title('Explicit Method: Error Plots for Stable $u_{h,\Delta t}(x,t)$');
xlb = xlabel('$x$');ylb = ylabel('$u(x)$');
set([t1 xlb ylb],'Interpreter','LaTeX');
xlb.FontSize=40;ylb.FontSize=40;t1.FontSize=30;
%print(gcf,'Diffusion_FTCS_errors.png','-dpng','-r500');



T = table(nx_vals',nt_vals',Error');
T.Properties.VariableNames = {'nx' 'nt' 'MaxError'}




animation_toggle = 1;
if ~animation_toggle
    return
end


% create figure
figure('name','diffusion u(x,t)','rend','painters','pos',[0 0 1000 1000]);
hold on
clf

%% animation
k=gca;
set(gca,'linewidth',3,'fontsize',20)
k.TickLength = [0.02 0.02];
j = 0;
% MyColor = flip([12 27 51;122 48 108;3 181 170;140 216 103;219 254 135]/255);
MyColor = flip([250 139 0;86 180 233;0 158 115;0 114 178;213 94 0]/255);

style = {'-';':';'--';'-.';'-'};
markers = {'+','o','s','x'};	
nt = 128;
[U,E,cond,t] = DiffusionFTCS(nx,nt,a,1,1,u_init,func_U);
for i = 1:4:54 % we start seeing instability for h = 1/32, dt = 1/128 around t = 45 
    j = j + 1;
    hold on
    k.YLimMode='manual';
    y = U(:,i);
    % axis([0 1 min(y)-5 max(y)+5])
    axis([0 1 -0.5 2])
    p(j) = plot(x_vec,y,'color','k','LineWidth',1); % unstable
    %p(j) = plot(x_vec,y,'color',MyColor(j,:),'LineWidth',4); % stable
    i
    %p(j) = plot(x_vec,y,'color',MyColor(j,:),'LineWidth',4,'Marker',markers{j},'LineStyle',style{j},'MarkerSize',8);
    pause(0.001);
    set(p(j),'Ydata',y);
    refreshdata
    drawnow
end

% u(x) profile

T = 1; % compare actual
xx = 0:0.001:1;
actual_soln = exp(-T)*sin(pi*xx)+exp(-4*T)*sin(2*pi*xx);
p2 = plot(xx,actual_soln,'color','k','linewidth',4,'LineStyle','-');%,'Marker','.','MarkerSize',15);
t2 = title(['Explicit Method: Diffusion ($u(x)$ profile) for $h = 1/$' num2str(nx) ' and $\Delta t = 1/$' num2str(nt) ' (' cond ')']);
xlb = xlabel('$x$');ylb = ylabel('$u(x)$');
set([t2 xlb ylb],'Interpreter','LaTeX');
xlb.FontSize=40;ylb.FontSize=40;t2.FontSize=25;

% lgd2 = legend([p p2],'$u(x,t_0)$',...
%     '$u(x,t_1)$',...
%     '$u(x,t_2$',...
%     '$u(x,t_3)$','$u(x,t_N)$');
% set(lgd2,'Interpreter','LaTeX');
% lgd2.FontSize=30;

print(gcf,'Diffusion_FTCS_ux_unstable.png','-dpng','-r500');


% theoretical U


% figure
% size(actual_u)
% size(u')
% [X,T] = meshgrid(x,t);
% surf(X,T,u')

% size(t)
% size(E)
% figure
% plot(t,E)

% plot(x,u(:,end),'ro','linewidth',1)



%% FTCS function
function [u,E,cond,t_vec] = DiffusionFTCS(nx,nt,a,xmax,tmax,u_init,func_U)

%function [u,E,cond] = MyDiffusionFTCS(nx,nt,a,xmax,tmax,u_init,func_U)
% 
% fprintf('nh = %f \n',nx) % number of spatial steps h
% fprintf('nt = %f \n',nt) % number of time steps
% a = pi^-2; % thermal diffusivity constant


x0 = 0;
x_vec = linspace(x0,xmax,nx);
h = xmax/(nx-1);

t0 = 0;
t_vec = linspace(t0,tmax,nt);
dt = tmax/(nt-1);


% Stability
if dt <= h^2/(2*a)
    cond = 'stable';
    %disp(['h = 1/' num2str(nx) ' and dt = 1/' num2str(nt) ' satisfy stability conditions'])
else
    cond = 'unstable';
    %disp(['h = 1/' num2str(nx) ' and dt = 1/' num2str(nt) ' DO NOT satisfy stability conditions'])
end

g0 = 0;
g1 = g0; % homogenous boundary conditions

% initialize u(x,t) array
u = zeros(nx,nt);
% ICs
u(1,1:end) = g0; % u(0,t) = g0(t)
u(nx,1:end) = g1; % u(1,t) = g1(t)
u(:,1) = u_init(x_vec);

for n = 1:nt-1
    for i = 2:nx-1
        u(i,n+1) = u(i,n) + (a*dt/h^2) * (u(i-1,n)-2*u(i,n)+u(i+1,n));
    end
end

u = u';
actual_u = func_U(x_vec,t_vec');
E = max(abs(u'-actual_u'));
u = u';
% u
% uf = u(:,end);
% figure
% plot(x,uf);
% keyboard

end



