clc
clear all
close all

u_init = @(x) sin(pi*x)+sin(2*pi*x);
func_U = @(x,t)(exp(-t)).*sin(pi*x) + (exp(-4*t)).*sin(2*pi*x);


%:imax % to collapse, cmd,+.  cmd,shft,+ to expand 
% diffusionCN(nt,nx,a,xmax,tmax,u_init,func_U)
nx = 17;
nt = 59;

a = pi^-2;
xmax = 1;
tmax = 1;
%% Theoretical Surface Plot
%     % Calculated U
%     figure(1)
%     subplot(1,imax,i)
%     hold on
%     colormap(jet)
%     AC = surf(X,T,U,'LineStyle','none');
%     view([0 90]); % birdseye view of surface plot
%     xlabel('x'),ylabel('t'),title('1D Diffusion using Crank-Nicholson Method')
%     colorbar
%% Crank Nicholson Surface Plot
%     figure('name','Diffusion - CN','rend','painters','pos',[0 0 900 900]);
%     clf
%     %subplot(3,imax+2,i) % 3 by 4 plot: 3 Explicit, 3 CN, 3 FEM, error plots of each
%     hold on
%     colormap(jet)
%     grid
%     surf(X,T,actual_U','LineStyle','none');
%         
%     view([-20,21])
% %     view([0 90]);
%     xlabel('x'),ylabel('t'),zlabel('u(x,t)'),title('Diffusion - Crank Nicholson Method')
%     set(gca,'linewidth',3,'fontsize',20)
%     colorbar
%print(gcf,'DiffusionSurf.png','-dpng','-r500');

%% Error Plot Comparison
tic;[u_FTCS,E_FTCS,cond] = MyDiffusionFTCS(nx,nt,a,xmax,tmax,u_init,func_U); t_FTCS = toc;
tic;[U,actual_U,E_CN,X,T] = MyDiffusionCN(nt,nx,a,xmax,tmax,u_init,func_U);t_CN = toc;
tic;[U_FEM,E_FEM] = MyDiffusionFEM(nt,nx,a,xmax,tmax,u_init,func_U); t_FEM = toc;
disp(cond)

% error plot
figure('name','Error Plots','rend','painters','pos',[0 0 900 900]);
clf
hold on
set(gca,'linewidth',3,'fontsize',20);
t = linspace(0,tmax,nt);

% FTCS stability check:
% delT = @(h) (h^2)/(2*pi^-2)
% floor(1/(0.9*delT(1/64))+1) % legal dt value
if strcmpi(cond,'stable')
    p_FTCS = plot(t,E_FTCS,'r--','LineWidth',3);
else
    p_FTCS = plot([0],[0],'r--','LineWidth',3);
end

p_CN = plot(t,E_CN,'b','LineWidth',3);
p_FEM = plot(t,E_FEM,'g:','LineWidth',3);
xlabel('t'),ylabel('max error'),title('Max u(x,t) error versus time (t)');
lgd = legend([p_FTCS,p_CN,p_FEM],'explicit','Crank-Nicholson','finite element');
%lgd = legend([p_CN,p_FEM],'Crank-Nicholson','finite element');

lgd.FontSize=30;
% print(gcf,'ErrorPlot.png','-dpng','-r500');

[t_FTCS t_CN t_FEM];

