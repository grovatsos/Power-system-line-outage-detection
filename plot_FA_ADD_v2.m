%8/31/2015
%Plots the ADD and FA for 3 bus and 14 bus systems


close all
clear all
clc


%%




%%
%14bus
FA=[1/24 1/4 1/2 1 2 7];
FA_cusum=FA.*24*3600*15;
A_cusum=log(FA_cusum./0.7265)/0.9545; %from excel best fit
ADD_cusum=0.0397.*A_cusum+0.0484; %from excel best fit
ADD_cusum=ADD_cusum./15;
%================================
% FA=[1/24 1/4 1/2 1 2 7];
FA_shewhart=FA.*24*3600*15;
A_shewhart=log(FA_shewhart./1.0962)/1.0682;

ADD_shewhart=1.0879*exp(0.026*A_shewhart);
ADD_shewhart=ADD_shewhart./15;
%================================
% FA=[1/24 1/4 1/2 1 2 7];
FA_meanshift=FA.*24*3600*15;
A_meanshift=log(FA_meanshift./0.5912)/1.6247;

ADD_meanshift=0.5931*exp(0.1042*A_meanshift);
ADD_meanshift=ADD_meanshift./15;
%================================


semilogx(FA,ADD_cusum, 'r');
set(gca,'XTick',FA);
set(gca,'XTickLabel',{});
set(gca,'XTickLabel',{'1/24','1/4','1/2','1','2','7'});
xl=xlabel('Mean Time to False Alarm [day]','fontsize',28,'fontname','times new roman');
x2=ylabel('Average Delay [sec]','fontsize',28,'fontname','times new roman');

hold on

semilogx(FA,ADD_shewhart, 'b');
set(gca,'XTick',FA);
set(gca,'XTickLabel',{});
set(gca,'XTickLabel',{'1/24','1/4','1/2','1','2','7'});
xl=xlabel('Mean Time to False Alarm [day]','fontsize',28,'fontname','times new roman');
x2=ylabel('Average Delay [sec]','fontsize',28,'fontname','times new roman');



semilogx(FA,ADD_meanshift, 'g');
set(gca,'XTick',FA);
set(gca,'XTickLabel',{});
set(gca,'XTickLabel',{'1/24','1/4','1/2','1','2','7'});
xl=xlabel('Mean Time to False Alarm [day]','fontsize',28,'fontname','times new roman');
x2=ylabel('Average Delay [sec]','fontsize',28,'fontname','times new roman');
legend('Cusum','Shewhart','Mean-Shift','location','northwest')
% ax = gca;
% ax.XTick = [1/24 1/4 1/2 1 2 7];
