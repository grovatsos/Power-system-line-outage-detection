close all
clear all
%generates the figures based on data from run_qcd_nonlinear code

figure('position', [20 80 1200 800]);
axes('position', [0.12 0.1 0.82 0.76]);

sample_per_sec = 15; %halved from 30 since we discard every other pt. due to iid issue

% loseline(75): (54,55)
load('L75out.mat');
PFI_g0
subplot(4,1,1);
set(gca, 'Position', [0.12,0.79,0.82,0.19]);
semilogx(beta,ADDvec/sample_per_sec, 'b--<', 'LineWidth', 2, 'MarkerSize', 13,...
        'MarkerFaceColor','b');
hold on; grid on;
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
set(gca,'XTick',beta);
set(gca,'XTickLabel',{});
set(gca,'fontsize',30,'fontname','times new roman');
ll=legend('Outage of (54,55)', 'Location','NorthWest');
set(ll,'Interpreter','latex');hold on
axis([0.7,250,3.51,5]);

yl=ylabel('Detection Delay [s]','fontsize',34,'fontname','times new roman');
set(yl,'Position',[0.4,1.6,1]);
set(yl,'Interpreter','latex');hold on

% loseline(91): (63,59)
load('L91out.mat');
PFI_g0
subplot(4,1,2);
set(gca, 'Position', [0.12,0.57,0.82,0.19]);
semilogx(beta,ADDvec/sample_per_sec, 'm--o', 'LineWidth', 2, 'MarkerSize', 13,...
        'MarkerFaceColor','m');
hold on; grid on;
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
set(gca,'XTick',beta);
set(gca,'XTickLabel',{});
set(gca,'fontsize',30,'fontname','times new roman');
ll=legend('Outage of (63,59)', 'Location','NorthWest');
set(ll,'Interpreter','latex');hold on
axis([0.7,250,0.121,0.16]);

% loseline(95): (64,65)
load('L95out.mat');
PFI_g0
subplot(4,1,3);
set(gca, 'Position', [0.12,0.35,0.82,0.19]);
semilogx(beta,ADDvec/sample_per_sec, 'r--s', 'LineWidth', 2, 'MarkerSize', 13,...
        'MarkerFaceColor','r');
hold on; grid on;
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
set(gca,'XTick',beta);
set(gca,'XTickLabel',{});
% set(gca,'YTick',[0.022 0.024]);
% set(gca,'YTickLabel',{'0.022','0.024'});
set(gca,'fontsize',30,'fontname','times new roman');
ll=legend('Outage of (64,65)', 'Location','NorthWest');
set(ll,'Interpreter','latex');hold on
axis([0.7,250,0.029,0.035]);
    
% loseline(102): (65,68)
load('L102out.mat');
PFI_g0
subplot(4,1,4);
set(gca, 'Position', [0.12,0.13,0.82,0.19]);
semilogx(beta,ADDvec/sample_per_sec, 'g--x', 'LineWidth', 2, 'MarkerSize', 13,...
        'MarkerFaceColor','r');
hold on; grid on;
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
set(gca,'XTick',beta);
set(gca,'XTickLabel',{'1/24','1/4','1/2','1','2','7'});
set(gca,'YTick',[0.0045 0.005 0.0055]);
set(gca,'YTickLabel',{'0.0045','0.005','0.0055'});
set(gca,'fontsize',30,'fontname','times new roman'); ...
   xl=xlabel('Mean Time to False Alarm [day]','fontsize',34,'fontname','times new roman');
   set(xl,'Interpreter','latex');
   ll=legend('Outage of (65,68)', 'Location','NorthWest');
   set(ll,'Interpreter','latex');hold on
axis([0.7,250,0.0045,0.0058]);
