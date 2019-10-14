% Plots Threshold, A vs. Average Detection Delay for line 74 outage of 118
% bus
clear all
close all
clc


A=[20 30 40 50 75 100 200 300 400 600 800];
ADD=[1.25 1.79 1.75 2.92 4.46 5.47 10.38 17.8 19.7 33.7 41.1];
%For False Alarm Plot
A=[7 10 14 18 25 40];
FA=1./[2.86 4.8 10.8 18 61 215];


figure('position', [20 80 1300 800]);
axes('position', [0.14 0.14 0.8 0.8]);
hold on; box on; grid on;
% title(['Double Line Outage ', num2str(linefault_dbl)])
p1=plot(A, FA,'--b', 'LineWidth', 3, 'Marker', '*', 'Markersize', 8);


set(gca,'fontsize',36,'fontname','times new roman'); ...
   xl=xlabel('Threshold $A$','fontsize',38,'fontname','times new roman');
   set(xl,'Interpreter','latex');
   yl=ylabel('Avg. False Alarm Rate','fontsize',38,'fontname','times new roman');
   set(yl,'Interpreter','latex');

%    set(l1,'FontSize',37, 'Interpreter','latex')
%    set(gca,'XTick',1.0487:0.0001:1.0492)
%    set(gca,'YTick',4952270:40:4952290)
%    axis([1.04901 1.0504 0.98 1.07])






