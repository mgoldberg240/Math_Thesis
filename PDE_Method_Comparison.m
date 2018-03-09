clc
clear all
u_init = @(x) sin(pi*x)+sin(2*pi*x);
func_U = @(x,t)(exp(-t)).*sin(pi*x) + (exp(-4*t)).*sin(2*pi*x);

imax = 3;
for i = 1:imax % to collapse, cmd,+.  cmd,shft,+ to expand 
    % diffusionCN(nt,nx,a,xmax,tmax,u_init,func_U)
    nt = 64*i;
    nx = 16*i;
    a = pi^-2;
    xmax = 1;
    tmax = 1;
    [U,actual_U,E,X,T] = MyDiffusionCN(nt,nx,a,xmax,tmax,u_init,func_U);
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
    figure(2)
    subplot(3,imax+2,i) % 3 by 4 plot: 3 Explicit, 3 CN, 3 FEM, error plots of each
    hold on
    colormap(jet)
    TC = surf(X,T,actual_U','LineStyle','none');
    view([0 90]);
    xlabel('x'),ylabel('t'),title('Theoretical U')
    colorbar

    % error plot
    figure(3)
    hold on
    t = linspace(0,tmax,nt);
    plot(t,E)
    xlabel('x'),ylabel('time'),title('Max U error versus time')
end
