%8/31/2015
%Plots the ADD and FA for 3 bus and 14 bus systems


close all
clear all
clc


%%
%3bus

FA=[1/24 1/4 1/2 1 2 7]; %days till false alarm
FA_cusum=FA.*24*3600*15; %number of samples till false alarm, assuming 15 samples per second (every other point for iid)
A_cusum=log(FA_cusum./4.812)/0.9323; %from excel best fit
ADD_cusum=0.2674.*A_cusum-0.602; %from excel best fit
ADD_cusum=ADD_cusum./15;

% FA=[1/24 1/4 1/2 1 2 7];
FA_shewhart=FA.*24*3600*15;
A_shewhart=log(FA_shewhart./3.5263)/1.0901;

ADD_shewhart=1.5756*exp(0.1293*A_shewhart);
ADD_shewhart=ADD_shewhart./15;

figure
title('3bus')
hold on
plot(FA,ADD_cusum, 'r');
plot(FA,ADD_shewhart, '*');

figure

plot(log(FA),ADD_cusum, 'r');

figure
plot(log(FA),log(ADD_shewhart), 'b');



%%
%14bus
FA=[1/24 1/4 1/2 1 2 7];
FA_cusum=FA.*24*3600*15;
A_cusum=log(FA_cusum./0.7265)/0.9545; %from excel best fit
ADD_cusum=0.0397.*A_cusum+0.0484; %from excel best fit
ADD_cusum=ADD_cusum./15;

% FA=[1/24 1/4 1/2 1 2 7];
FA_shewhart=FA.*24*3600*15;
A_shewhart=log(FA_shewhart./1.0962)/1.0682;

ADD_shewhart=1.0879*exp(0.026*A_shewhart);
ADD_shewhart=ADD_shewhart./15;

figure
hold on
title('14bus')
plot(FA,ADD_cusum, 'r');
plot(FA,ADD_shewhart, '*');

figure

semilogx(FA,ADD_cusum, 'r');
set(gca,'XTick',FA);
set(gca,'XTickLabel',{});
set(gca,'XTickLabel',{'1/24','1/4','1/2','1','2','7'});
xl=xlabel('Mean Time to False Alarm [day]','fontsize',34,'fontname','times new roman');

figure
loglog(FA,ADD_shewhart, '*');
set(gca,'XTick',FA);
set(gca,'XTickLabel',{});
set(gca,'XTickLabel',{'1/24','1/4','1/2','1','2','7'});
xl=xlabel('Mean Time to False Alarm [day]','fontsize',34,'fontname','times new roman');
% ax = gca;
% ax.XTick = [1/24 1/4 1/2 1 2 7];
