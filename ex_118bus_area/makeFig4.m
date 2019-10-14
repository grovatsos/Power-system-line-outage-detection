close all
clear all

figure('position', [20 80 1200 800]);
axes('position', [0.12 0.1 0.82 0.76]);

% loseline(12): (3,12)
load('L12out.mat');
PFI_g0
subplot(4,1,1);
set(gca, 'Position', [0.12,0.79,0.82,0.19]);
semilogx(beta,ADDvec/30, 'b--<', 'LineWidth', 2, 'MarkerSize', 13,...
        'MarkerFaceColor','b');
hold on;
% semilogx(, 'b--<', 'LineWidth', 2, 'MarkerSize', 13,...
%         'MarkerFaceColor','b');
grid on;
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
set(gca,'XTick',beta);
set(gca,'XTickLabel',{});
set(gca,'fontsize',30,'fontname','times new roman');
ll=legend('Outage of (3,12)', 'Location','NorthWest');
set(ll,'Interpreter','latex');hold on
axis([0.7,250,0.45,0.7]);

yl=ylabel('Detection Delay [s]','fontsize',34,'fontname','times new roman');
set(yl,'Position',[0.4,0.12,1]);
set(yl,'Interpreter','latex');hold on


% loseline(30): (25,27)
load('L30out.mat');
PFI_g0
subplot(4,1,2);
set(gca, 'Position', [0.12,0.57,0.82,0.19]);
semilogx(beta,ADDvec/30, 'm--o', 'LineWidth', 2, 'MarkerSize', 13,...
        'MarkerFaceColor','m');
hold on; grid on;
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
set(gca,'XTick',beta);
set(gca,'XTickLabel',{});
set(gca,'fontsize',30,'fontname','times new roman');
ll=legend('Outage of (25,27)', 'Location','NorthWest');
set(ll,'Interpreter','latex');hold on
axis([0.7,250,5.5/30,8/30]);


% loseline(14): (11,13)
load('L14out.mat');
PFI_g0
subplot(4,1,3);
set(gca, 'Position', [0.12,0.35,0.82,0.19]);
semilogx(beta,ADDvec/30, 'r--s', 'LineWidth', 2, 'MarkerSize', 13,...
        'MarkerFaceColor','r');
hold on; grid on;
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
set(gca,'XTick',beta);
set(gca,'XTickLabel',{});
set(gca,'YTick',[0.022 0.024]);
set(gca,'YTickLabel',{'0.022','0.024'});
set(gca,'fontsize',30,'fontname','times new roman');
ll=legend('Outage of (11,13)', 'Location','NorthWest');
set(ll,'Interpreter','latex');hold on
axis([0.7,250,0.021,0.026]);
    
% loseline(7): (8,5)
load('L7out.mat');
PFI_g0
subplot(4,1,4);
set(gca, 'Position', [0.12,0.13,0.82,0.19]);
semilogx(beta,ADDvec/30, 'g--x', 'LineWidth', 2, 'MarkerSize', 13,...
        'MarkerFaceColor','r');
hold on; grid on;
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
set(gca,'XTick',beta);
set(gca,'XTickLabel',{'1/24','1/4','1/2','1','2','7'});
set(gca,'YTick',[0.0025 0.003]);
set(gca,'YTickLabel',{'0.0025','0.003'});
set(gca,'fontsize',30,'fontname','times new roman'); ...
   xl=xlabel('Mean Time to False Alarm [day]','fontsize',34,'fontname','times new roman');
   set(xl,'Interpreter','latex');
   ll=legend('Outage of (8,5)', 'Location','NorthWest');
   set(ll,'Interpreter','latex');hold on
axis([0.7,250,0.075/30,0.099/30]);
