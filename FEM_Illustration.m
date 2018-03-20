clc
clear
close all

% sketch function in question
L = 1;
x = linspace(0,L,1000);
F = sin(pi*x)+sin(2*pi*x);


% create unit tents
N = 10; % number of tents
xTents = linspace(0,L,N); xstep = xTents(1);
FTents = sin(pi*xTents)+sin(2*pi*xTents);


yTents = ones(length(xTents),1);
yTents(2:2:end) = 0;
yTentsP = abs(yTents-1)';

figure('name','Finite Element Approximation of u(x)','rend','painters','pos',[0 0 900 900]);
clf
hold on
% set plot limits
lims = num2cell([-0.1,1.1,-0.5,2.1]);
[xmin,xmax,ymin,ymax] = deal(lims{:}); % assigned four variables in two lines
axis([xmin xmax ymin ymax])

% plot F and tents
Fplot = plot(x,F,'blue','Linewidth',4,'DisplayName','u(x) = sin(pi*x)+sin(2*pi*x)'); % plot F(x)
Tplot = plot(xTents,yTents,'k','Linewidth',3,'DisplayName','Tent Functions'); % plot some tents
plot(xTents,yTentsP,'k','Linewidth',3) % plot the remaining tents

% plot vertical lines from tents to F
for x = 1:length(xTents)
    plot([xTents(x) xTents(x)],[0 FTents(x)],'g','Linewidth',2) % plot green vertical lines from F to 0
end

% plot FEM approximation of F
FEM = plot(xTents',FTents','r:','Linewidth',4,'DisplayName','Finite Element Approximation of u(x)');

% plot x and y axis
% plot([0 1],[0 0],'k','LineWidth',2)
% plot([0 0],[0 2],'k','LineWidth',2)

% label figure
xlabel('x');ylabel('u(x)')
title('Finite Element Approximation of u(x)')
lgd = legend([Fplot Tplot FEM],{'u(x) = sin(\pix)+sin(2\pix)','Tent Functions','Finite Element Approximation of u(x)'});
lgd.FontSize=20;
set(gca,'linewidth',3,'fontsize',20)

% tr = [0,1,0];
% draw_tr = area([xstep,xstep+1/(N+1),xstep+2/(N+1)],tr);
% draw_tr.FaceAlpha = 0.5;
% draw_tr.FaceColor = [253 208 35]/255;

% plot text "phi"
% dim = [0.26 0.3 0.15 0.16];
% str = {'$\phi_1$'};
% annotation('textbox',dim,...
%     'String',str,...
%     'Interpreter','Latex',...
%     'FitBoxToText','on',...
%     'Linestyle','None',...
%     'Fontsize',48,...
%     'FontWeight','bold');

print(gcf,'FEM_Illustration.png','-dpng','-r500');
