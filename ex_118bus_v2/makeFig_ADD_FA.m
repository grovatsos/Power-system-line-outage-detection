close all
clear all
%generates the figures based on data from run_qcd_nonlinear code

figure('position', [20 80 1200 800]);
axes('position', [0.12 0.1 0.82 0.76]);

sample_per_sec = 15;
sample_per_day = 15*60*60*24;

A=[5:5:80];
ADD_Cusum=0.0525*A+0.267;
FA_Cusum=1.479*exp(0.1318*A);


ADD_Shrewhart=1.3189*exp(0.0314*A);
FA_Shrewhart=0.4689*exp(0.3545*A);

% set(gca, 'Position', [0.12,0.79,0.82,0.19]);
subplot(2,1,1);
semilogx(FA_Cusum/sample_per_sec,ADD_Cusum/sample_per_sec, 'b--*', 'LineWidth', 2, 'MarkerSize', 7,...
        'MarkerFaceColor','b');
hold on; grid on;
   
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
% set(gca,'XTick',FA_Cusum);
% set(gca,'XTickLabel',{});
set(gca,'fontsize',24,'fontname','times new roman');
ll=legend('CuSum: Line 74 Outage', 'Location','NorthWest');
set(ll,'Interpreter','latex');
% axis([0.7,250,0,0.4]);

% yl=ylabel('Detection Delay [s]','fontsize',33,'fontname','times new roman');
% ylabel_handle=get(gca,'YLabel');
% set(ylabel_handle,'Position',get(ylabel_handle,'Position') - [0 .4 0]);

% set(yl,'Position',[0.9,4.6,1]);

subplot(2,1,2);
semilogx(FA_Shrewhart/sample_per_sec,ADD_Shrewhart/sample_per_sec, 'r--*', 'LineWidth', 2, 'MarkerSize', 7,...
        'MarkerFaceColor','r');
hold on; grid on;
set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
% set(gca,'XTick',FA_Shrewhart);
% set(gca,'XTickLabel',{});
set(gca,'fontsize',24,'fontname','times new roman');
ll=legend('Shewhart: Line 74 Outage', 'Location','NorthWest');
set(ll,'Interpreter','latex');
% axis([0.7,250,3.51,5]);




% set(gca,'fontsize',30,'fontname','times new roman');
xl=xlabel('Mean Time to False Alarm [s]','fontsize',33,'fontname','times new roman');
set(xl,'Interpreter','latex');
%    ll=legend('Outage of (65,68)', 'Location','NorthWest');
%    set(ll,'Interpreter','latex');hold on
axis([0.01,10^12,0,1.5]);
yl=ylabel('Detection Delay [s]','fontsize',33,'fontname','times new roman');
set(yl,'Interpreter','latex');
ylabel_handle=get(gca,'YLabel');
set(ylabel_handle,'Position',get(ylabel_handle,'Position') + [0 1 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_max=log(FA_Cusum(end)/0.4689)/0.3545;
A_min=log(FA_Cusum(1)/0.4689)/0.3545;
A_Shrewhart=linspace(A_min,A_max,numel(A))

ADD_Shrewhart=1.3189*exp(0.0314*A_Shrewhart)
FA_Shrewhart=0.4689*exp(0.3545*A_Shrewhart)

figure('position', [20 80 1200 800]);
axes('position', [0.12 0.1 0.82 0.76]);
% set(gca, 'Position', [0.12,0.79,0.82,0.19]);
s1=semilogx(FA_Cusum/sample_per_sec, ADD_Cusum/sample_per_sec, 'b--*', 'LineWidth', 2, 'MarkerSize', 7,...
        'MarkerFaceColor','b');
hold on; grid on;
s2=semilogx(FA_Shrewhart/sample_per_sec,ADD_Shrewhart/sample_per_sec, 'r--*', 'LineWidth', 2, 'MarkerSize', 7,...
        'MarkerFaceColor','r');

set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','XMinorTick','off');
% set(gca,'XTick',FA_Cusum);
% set(gca,'XTickLabel',{});
set(gca,'fontsize',24,'fontname','times new roman');
ll=legend([s1,s2], 'CuSum', 'Shewhart', 'Location','NorthWest');
set(ll,'Interpreter','latex');
% axis([0.7,250,0,0.4]);

% yl=ylabel('Detection Delay [s]','fontsize',33,'fontname','times new roman');
% ylabel_handle=get(gca,'YLabel');
% set(ylabel_handle,'Position',get(ylabel_handle,'Position') - [0 .4 0]);

% set(yl,'Position',[0.9,4.6,1]);


xl=xlabel('Mean Time to False Alarm [s]','fontsize',33,'fontname','times new roman');
set(xl,'Interpreter','latex');
axis([0.1,10^4,0,0.5]);
yl=ylabel('Detection Delay [s]','fontsize',33,'fontname','times new roman');
set(yl,'Interpreter','latex');
% ylabel_handle=get(gca,'YLabel');
% set(ylabel_handle,'Position',get(ylabel_handle,'Position') + [0 1 0]);

