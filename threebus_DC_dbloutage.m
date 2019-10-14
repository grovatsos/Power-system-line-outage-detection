%%%fourteenbus_DC is the one that runs the 14bus example, generates sample paths of Wn from CuSum algorithm 
%7/16/2014: Automate some of the code so it can be applied to larger
%systems
%7/23/2014: run sample paths for simulations, see second block of this code
%8/5/2014: Note, going from 3 bus system to 14 bus system, the variance on P and Q has been changed,
           %Some lines are better for simulating outages as they affect theta significantly
%10/28/2014: This code checks ADD by setting individual thresholds for each line outage scaled according to their KL divergence
%11/18/2014: Adapted from 3bus to 14bus system. 
%3/17/2015: Adapted to include double line outages
%4/25/2015: Code used to test Matpower pf. 3 bus system can't have double outage
% This code allows for PMUs at a subset of buses in addition to computing sample paths

clear all
close all
% clc

[bus, line] = buslinedata(3);
nbus = size(bus,1);
nline = size(line,1);
%%%%%%%%%%%
% line(29,2) = 22; %This is only for the 30 bus system to match diagram
%%%%%%%%%%%
%%generate a random covariance matrix with only diagonal entries

[soln, ybus, J] = get_pf_v4(bus,line);
% V0 = soln(1:end,2);
theta0 = soln(1:end,3).*pi/180; %in radians
J0 = zeros(nbus,nbus);
for i = 1:nline
    from = line(i,1);
    to = line(i,2);
    Xl = line(i,4); 
    J0(from,to) = -1/Xl;
    J0(to,from) = -1/Xl;
    J0(from,from) = J0(from,from) + 1/Xl;
    J0(to,to) = J0(to,to) + 1/Xl;
end

% %"internal system" A1
% J11 = J0(A1,A1);    
% J11(:,1) = []; %remove slack bus equations
% J11(1,:) = []; %bus 1 is assumed slack bus
% J12 = J0(A1,A2);    
% J12(1,:) = [];

%entire system
J0_orig=J0;
J0(:,1) = [];   J0(1,:) = []; %remove slack bus

%Compute single outage KL divergences

%use 0.05 for 14 bus, 0.1 for 3 bus, this does affect KL divergence so in theory, it shouldn't affect ADD
%but, if it hits some limit in power flow and goes to a different operating point the H matrices may change completely 
Psig = 0.3; %0.3 for 3bus system, 0.03 for 14 bus system
Qsig=0.02;

B0 = inv(J0);
sigma0 = B0*(2*Psig^2*eye(nbus-1))*B0'; %M*Sigma*M' in paper, this is all A1 subsystem

loseline_single = 1:nline;
B_single = cell(numel(loseline_single),1);
sigma_single = cell(numel(loseline_single),1);
% mu = cell(numel(loseline),1);
% klmean_single = zeros(numel(loseline),1);
KL_single = zeros(numel(loseline_single),2); %KL divergence and associated line loss

%Compute KL divergence for losing a line
for i = 1:numel(loseline_single)
    from = line(loseline_single(i),1);
    to = line(loseline_single(i),2);
    xij = line(loseline_single(i),4);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
    fij(1) = []; %remove slack bus
    dJ = -1/xij*(fij*fij');
    Ji = J0 + dJ; %H matrix in paper
    
    if rank(Ji) < size((Ji),1)
        disp(['critical single line ' num2str(loseline_single(i))]);
        disp(['line ' num2str(loseline_single(i)), ' has buses ' num2str(from), ' and ' num2str(to), sprintf('\n')]);
        B_single{i}=NaN;
        sigma_single{i}=NaN;
        KL_single(i,1) = NaN; 
        KL_single(i,2) = loseline_single(i); 

        continue
    end
    
    B_single{i} = inv(Ji); %Compute each M corresponding to each line outage in A1;
    sigma_single{i} = Ji\(2*Psig^2*eye(nbus-1))/(Ji'); %make it symmetric to avoid numerical rounding issues
    aa = triu(sigma_single{i});
    sigma_single{i} = aa + triu(aa,1)';

    KL_single(i,1) = 1/2 *(trace(sigma0\sigma_single{i}) - ...
        (length(sigma0)) + log(det(sigma0/sigma_single{i}))); 
    KL_single(i,2) = loseline_single(i); 
end

   
%Lowest KL divergence dictates the misdetection error. We want KL divergences to be big
KL_single_sorted=sortrows(KL_single,1); %sort KL divergnce from min to max as well as their associated lines
%disp(['min KL divergence for Area 1 single outage = ' num2str(KL_single_sorted(1,1))]);
critical_lines_single=find(isnan(KL_single)); %lines that if lost, creates islanding situations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Limited PMU case
remove_PMU_bus=[]; %buses that do not have PMUs, bus 1 is slack, so remove PMU can only remove from 2 or 3

% Compute the KL divergences with certain bus' PMUs removed
PMU_bus=setdiff(1:nbus, [1 remove_PMU_bus]); %bus that has the PMU
B0(remove_PMU_bus-1,:)=[];
sigma0 = B0*(2*Psig^2*eye(nbus-1))*B0';
KL_single_limitedPMUs = zeros(numel(loseline_single),2); %KL divergence and associated line loss with less than full set of PMUs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute the new KL divergences with limited number of PMUs for single line outages
for i = 1:numel(loseline_single)
    if ismember(i, critical_lines_single)
         KL_single_limitedPMUs(i,1) = NaN;
         KL_single_limitedPMUs(i,2)=loseline_single(i);
        continue
    end
    Bi=B_single{i};
    Bi(remove_PMU_bus-1,:)=[];
    B_single{i} = Bi;
    sigma_single{i} = Bi*(2*Psig^2*eye(nbus-1))*Bi'; %make it symmetric to avoid numerical rounding issues
    aa = triu(sigma_single{i});
    sigma_single{i} = aa + triu(aa,1)';
    KL_single_limitedPMUs(i,1) = KL_compute(sigma_single{i}, sigma0); 
    KL_single_limitedPMUs(i,2)=loseline_single(i);
end


KL_single_limitedPMUs_sorted=sortrows(KL_single_limitedPMUs,1); %sort KL divergnce from min to max as well as their associated lines


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create sample paths for a single line fault

% seed=4;
% rng(seed)
% linefault=datasample(loseline,1); %choose a random line to fault
linefault_single=3;
faulttime=1;%fault time
if (ismember(linefault_single, critical_lines_single)|| linefault_single<1 || linefault_single > numel(loseline_single))
    error('choose different line to fault')
end
num_samplepaths=1;
y=theta0;
ADD=Inf;
%run sample paths
%scale thresholds according to KL divergences
A_single=(100*KL_single_limitedPMUs(:,1)/max(KL_single_limitedPMUs(:,1)))';
A_single=100*ones(1, numel(loseline_single)); %Threshold

mu0=zeros(1, numel(theta0));
theta0_fault=cell(numel(loseline_single,1));
for i=1:numel(loseline_single)
    if ismember(i, critical_lines_single)
        theta0_fault{i}=NaN;
        continue
    end
    line_temp=line;
    line_temp(i,:)=[]; %remove faulted line
    [soln_fault, ybus_fault, J_fault] = get_pf_v4(bus,line_temp);
    theta0_fault{i} = soln_fault(1:end,3).*pi/180;
end

for pathId=1:num_samplepaths 
    tao_single=1; %time
    f_stop=0;
    f_single=[]; %Specifies which single line outage streams have crossed threshold 
    f_dbl=[];
    dyvec = [];
    Wn_single=zeros(1, numel(loseline_single)); %CuSum, one for each line single line outage 
    Wnvec_single = [];

    
    %initialize first instant of y. If the simulation starts off with fault, then ynom is faulted y value. 
    if tao_single >= faulttime 
       y=theta0_fault{linefault_single};
    end    
    
    while f_stop==0
%       rng(seed)
        P=bus(2:end,6); %Negative signs are to match earlier version of code. Note that Pload and Qload may become negative meaning bus is generating power. 
        Q=bus(2:end,7);
        P_sample=P+randn(numel(P),1)*Psig;
        Q_sample=Q+randn(numel(Q),1)*Qsig;
         
        if tao_single < faulttime %no fault
            bus_sample=bus;
            bus_sample(2:end,6)=P_sample;
            bus_sample(2:end,7)=Q_sample;
            [soln_sample, ybus_sample, J_sample] = get_pf_v4(bus_sample,line);
%             V0 = soln(1:end,2);
            theta_sample = soln_sample(1:end,3).*pi/180; %in radians
            ynew=theta_sample;
            dy = ynew'-y';
            dy = dy(2:end);
            y = ynew;
        else
            bus_sample=bus;
            bus_sample(2:end,6)=P_sample;
            bus_sample(2:end,7)=Q_sample;
            line_temp=line;
            line_temp(linefault_single,:)=[];
            
            [soln_sample, ybus_sample, J_sample] = get_pf_v4(bus_sample,line_temp);
%             V0 = soln(1:end,2);
            theta_sample = soln_sample(1:end,3).*pi/180; %in radians
            ynew=theta_sample;
            dy = ynew'-y';
            dy = dy(2:end);
            y = ynew;
           
        end
        dyvec = [dyvec; dy];

        %Compute the CUSUM statistics for each possible post-change scenario
        for i=1:numel(loseline_single)
            if ismember(i, critical_lines_single)
                continue
            end
            Wn_single(i) = subplus(Wn_single(i) + log(mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma_single{i})/mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma0)));
        end
        Wnvec_single = [Wnvec_single; Wn_single];

            
        
        %Check if any one of them is above A
        if sum(subplus(Wn_single-A_single))>0
            [temp_single max_single]=max(Wn_single-A_single);
            Sorted_single=sort(subplus(Wn_single-A_single), 'descend');
%               find(temp(1)==(Wn_single-A_single))
            f_single=find(subplus(Wn_single-A_single)); 
            disp(['single line outage streams ', num2str(f_single), ' crossed the treshold ', num2str(A_single(f_single))]);       
            if pathId==1 
                ADD = tao_single;
            else
                ADD = ((pathId-1)*ADD + tao_single)/pathId; %running average to compute ADD over many sample paths
            end
            pathId
            f_stop==1;
            break; %needed so tao_single doesn't increase  
        end 
        %Update the stopping variable till we actually stop
        tao_single=tao_single+1;
        
        if tao_single>10000
            error('timed out')
        end
    end
end
ADD=ADD-faulttime+1;
disp(['single line ' num2str(linefault_single), ' fault ADD=' num2str(ADD)]);

figure('position', [20 80 1300 800]);
axes('position', [0.14 0.14 0.8 0.8]);
hold on; box on; grid on;
title(['Single Line Outage ', num2str(linefault_single)])
p1=plot(1:tao_single, Wnvec_single(1:tao_single,1)','--b', 'LineWidth', 2, 'Marker', '*', 'Markersize', 6);
p2=plot(1:tao_single, Wnvec_single(1:tao_single,2)','--g', 'LineWidth', 2,'Marker', 'o','Markersize', 6);
p3=plot(1:tao_single, Wnvec_single(1:tao_single,3)','--r', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% p4=plot(1:tao_single, Wnvec_single(1:tao_single,6)','--c', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% p5=plot(1:tao_single, Wnvec_dbl(1:tao_single,max_dbl)','--m', 'LineWidth', 2,'Marker', 's','Markersize', 6);
l1=legend([p1, p2, p3], '$W_{1}$','$W_{2}$','$W_{3}$');

set(gca,'fontsize',36,'fontname','times new roman'); ...
   xl=xlabel('$k$','fontsize',38,'fontname','times new roman');
   set(xl,'Interpreter','latex');
   yl=ylabel('$W[k]$','fontsize',38,'fontname','times new roman');
   set(yl,'Interpreter','latex');
   set(l1,'Interpreter','latex');
   set(l1,'FontSize',37, 'Interpreter','latex')
%    set(gca,'XTick',1.0487:0.0001:1.0492)
%    set(gca,'YTick',4952270:40:4952290)
%    axis([1.04901 1.0504 0.98 1.07])

% %Plot all streams 
% figure('position', [20 80 1300 800]);
% axes('position', [0.14 0.14 0.8 0.8]);
% hold on; box on; grid on;
% title(['All streams for Single Line Outage ', num2str(linefault_single)])
% plot(1:tao_single, Wnvec_single(1:tao_single, :)', '--');
% plot(1:tao_single, Wnvec_dbl(1:tao_single, :)','-');

% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Create sample paths for double line outages
% % seed=4;
% linefault_dbl=[6 9]; %choose the two lines that are outaged [3 19] [3 20] show good paths, [6 9] produce example of case where line 5 single outage rate is faster than line 2 
% index=find((logical(loseline_dbl(:,1)==linefault_dbl(1))&logical(loseline_dbl(:,2)==linefault_dbl(2)))|(logical(loseline_dbl(:,1)==linefault_dbl(2))&logical(loseline_dbl(:,2)==linefault_dbl(1))));
% faulttime=1;%fault time
% if ismember(linefault_dbl, critical_lines_dbl, 'rows')
%     error('choose different two lines to fault')
% end
% num_samplepaths=1;
% y=theta0;
% ADD=Inf;
% %run sample paths
% %scale thresholds according to KL divergences
% A_single=(100*KL_single_limitedPMUs(:,1)/max(KL_single_limitedPMUs(:,1)))';
% A_single=300*ones(1, numel(loseline_single)); %Threshold
% A_dbl=(100*KL_dbl_limitedPMUs(:,1)/max(KL_dbl_limitedPMUs(:,1)))';
% A_dbl=300*ones(1, length(loseline_dbl)); %Threshold
% 
% mu0=zeros(1, numel(theta0));
% theta0_fault=cell(length(loseline_dbl), 1);
% for i=1:length(loseline_dbl)
%     if ismember(loseline_dbl(i,:), critical_lines_dbl, 'rows')
%         theta0_fault{i}=NaN;
%         continue
%     end
%     line_temp=line;
%     line_temp(linefault_dbl,:)=[]; %remove faulted line
%     [soln_fault, ybus_fault, J_fault] = get_pf_v4(bus,line_temp);
%     theta0_fault{i} = soln_fault(1:end,3).*pi/180;
% end
% 
% for pathId=1:num_samplepaths 
%     tao_dbl=1; %time
%     f_stop=0;
%     f_single=[]; %Specifies which single line outage streams have crossed threshold 
%     f_dbl=[];
%     dyvec = [];
%     Wn_single=zeros(1, numel(loseline_single)); %CuSum, one for each line single line outage 
%     Wnvec_single = [];
%     
%     Wn_dbl=zeros(1, length(loseline_dbl));
%     Wnvec_dbl = [];
%     
%     %initialize first instant of y. If the simulation starts off with fault, then ynom is faulted y value. 
%     if tao_dbl >= faulttime 
%       y=theta0_fault{index};
%     end    
%     
%     while f_stop==0
% %       rng(seed)
%         P=-bus(2:end,6); %Negative signs are to match earlier version of code. Note that Pload and Qload may become negative meaning bus is generating power. 
%         Q=-bus(2:end,7);
%         P_sample=P+randn(numel(P),1)*Psig;
%         Q_sample=Q+randn(numel(Q),1)*Qsig;
%          
%         if tao_dbl < faulttime %no fault
%             bus_sample=bus;
%             bus_sample(2:end,6)=-P_sample;
%             bus_sample(2:end,7)=-Q_sample;
%             [soln_sample, ybus_sample, J_sample] = get_pf_v4(bus_sample,line);
% %             V0 = soln(1:end,2);
%             theta_sample = soln_sample(1:end,3).*pi/180; %in radians
%             ynew=theta_sample;
%             dy = ynew'-y';
%             dy = dy(2:end);
%             y = ynew;
%         else
%             bus_sample=bus;
%             bus_sample(2:end,6)=-P_sample;
%             bus_sample(2:end,7)=-Q_sample;
%             line_temp=line;
%             line_temp(linefault_dbl,:)=[];
%             [soln_sample, ybus_sample, J_sample] = get_pf_v4(bus_sample,line_temp);
% %             V0 = soln(1:end,2);
%             theta_sample = soln_sample(1:end,3).*pi/180; %in radians
%             ynew=theta_sample;
%             dy = ynew'-y';
%             dy = dy(2:end);
%             y = ynew;
%            
%         end
%         dyvec = [dyvec; dy];
% 
%         %Compute the CUSUM statistics for each possible post-change scenario
%         for i=1:numel(loseline_single)
%             if ismember(i, critical_lines_single)
%                 continue
%             end
%             Wn_single(i) = subplus(Wn_single(i) + log(mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma_single{i})/mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma0)));
%         end
%         Wnvec_single = [Wnvec_single; Wn_single];
%         
%         for i = 1:length(loseline_dbl)   
%             if ismember(loseline_dbl(i,:), critical_lines_dbl, 'rows') %is part of critical double line outage
%                 continue
%             end
%             Wn_dbl(i) = subplus(Wn_dbl(i) + log(mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma_dbl{loseline_dbl(i,1), loseline_dbl(i,2)})/mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma0)));
%         end
%         Wnvec_dbl = [Wnvec_dbl; Wn_dbl];
%             
%         
%         %Check if any one of them is above A
%         if (sum(subplus(Wn_single-A_single))>0) ||(sum(subplus(Wn_dbl-A_dbl))>0)
%             [temp_single max_single]=max(Wn_single-A_single);
%             [temp_dbl max_dbl]=max(Wn_dbl-A_dbl);
%             if sum(subplus(Wn_single-A_single))>0
%                 Sorted_single=sort(subplus(Wn_single-A_single), 'descend');
% %               find(temp(1)==(Wn_single-A_single))
%                 f_single=find(subplus(Wn_single-A_single)); 
%                 disp(['single line outage streams ', num2str(f_single), ' crossed the treshold ', num2str(A_single(f_single))]);
%                 warning('misdetection of double line fault')
%             end
%             if sum(subplus(Wn_dbl-A_dbl))>0
%                 Sorted_dbl=sort(subplus(Wn_dbl-A_dbl), 'descend');
% %               find(temp(1)==(Wn_dbl-A_dbl))
%                 f_dbl=find(subplus(Wn_dbl-A_dbl)); 
%                 disp('these dbl line outage streams crossed the treshold'); loseline_dbl(f_dbl, 1:2)' 
%                 disp('the corresponding thresholds are'); A_dbl(f_dbl)
%             end
%             if pathId==1 
%                 ADD = tao_dbl;
%             else
%                 ADD = ((pathId-1)*ADD + tao_dbl)/pathId; %running average to compute ADD over many sample paths
%             end
%             pathId
%             f_stop==1;
%             break; %needed so tao_single doesn't increase  
%         end 
%         %Update the stopping variable till we actually stop
%         tao_dbl=tao_dbl+1;
%         
%         if tao_dbl>1000
%             error('timed out')
%         end
%     end
% end
% ADD=ADD-faulttime+1;
% disp(['double line ' num2str(linefault_dbl), ' fault ADD=' num2str(ADD)]);
% 
% figure('position', [20 80 1300 800]);
% axes('position', [0.14 0.14 0.8 0.8]);
% hold on; box on; grid on;
% % title(['Double Line Outage ', num2str(linefault_dbl)])
% p1=plot(1:tao_dbl, Wnvec_dbl(1:tao_dbl, index )','--b', 'LineWidth', 2, 'Marker', '*', 'Markersize', 6);
% p2=plot(1:tao_dbl, Wnvec_dbl(1:tao_dbl, 3)','--g', 'LineWidth', 2,'Marker', 'o','Markersize', 6);
% p3=plot(1:tao_dbl, Wnvec_dbl(1:tao_dbl, 4)','--c', 'LineWidth', 2,'Marker', 'o','Markersize', 6);
% p4=plot(1:tao_dbl, Wnvec_single(1:tao_dbl, linefault_dbl(1))','-r', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% p5=plot(1:tao_dbl, Wnvec_single(1:tao_dbl, linefault_dbl(2))','-m', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% p6=plot(1:tao_dbl, Wnvec_single(1:tao_dbl, 3)','-k', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% % p7=plot(1:tao_dbl, Wnvec_single(1:tao_dbl, 12)','-k', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% % p5=plot(1:tao_dbl, Wnvec_single(1:tao_dbl, max_single)','--k', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% l1=legend([p1, p2, p3, p4, p5, p6], '$W_{(dbl)}$','$W_{(dbl 3)}$', '$W_{(dbl 4)}$','$W_{(6)}$','$W_{(9)}$', '$W_{(3)}$');
% 
% 
% set(gca,'fontsize',36,'fontname','times new roman'); ...
%    xl=xlabel('$k$','fontsize',38,'fontname','times new roman');
%    set(xl,'Interpreter','latex');
%    yl=ylabel('$W[k]$','fontsize',38,'fontname','times new roman');
%    set(yl,'Interpreter','latex');
%    set(l1,'Interpreter','latex');
%    set(l1,'FontSize',37, 'Interpreter','latex')
% %    set(gca,'XTick',1.0487:0.0001:1.0492)
% %    set(gca,'YTick',4952270:40:4952290)
% %    axis([1.04901 1.0504 0.98 1.07])
% 
% % %Plot all streams 
% % figure('position', [20 80 1300 800]);
% % axes('position', [0.14 0.14 0.8 0.8]);
% % hold on; box on; grid on;
% % title(['All streams for Double Line Outage ', num2str(linefault_dbl)])
% % plot(1:tao_dbl, Wnvec_single(1:tao_dbl, :)', '--');
% % plot(1:tao_dbl, Wnvec_dbl(1:tao_dbl, :)','-');
%    
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Compute KL divergences for Single and Double Line outages and also some False Isolation parameters
% 
% %single line outage, single line outage streams
% E_single_single=zeros(numel(loseline_single), numel(loseline_single)); %E_single_single(i,j) gives expected value of log(fj/f0) given line i went out 
% for i=1:numel(loseline_single)
%     if ismember(loseline_single(i), critical_lines_single)
%         E_single_single(i,:)=NaN;
%         continue
%     end
%     for j=1:numel(loseline_single)
%         if ismember(j, critical_lines_single)
%             E_single_single(i,j)=NaN;
%             continue
%         end
%         E_single_single(i,j)=KL_compute(sigma_single{i}, sigma0)-KL_compute(sigma_single{i}, sigma_single{j}); 
%     end
% end
% 
% %single line outage, double line outage streams
% E_single_dbl=zeros(numel(loseline_single), length(loseline_dbl)); %E_single_dbl(i,j) gives expected value of log(fj/f0) given line i went out, fj is double outage
% for i=1:numel(loseline_single)
%     if ismember(loseline_single(i), critical_lines_single)
%         E_single_dbl(i,:)=NaN;
%         continue
%     end
%     for j=1:length(loseline_dbl)
%         if ismember(loseline_dbl(j,:), critical_lines_dbl, 'rows')
%             E_single_dbl(i,j)=NaN;
%             continue
%         end
%         E_single_dbl(i,j)=KL_compute(sigma_single{i}, sigma0)-KL_compute(sigma_single{i}, sigma_dbl{loseline_dbl(j,1),loseline_dbl(j,2)}); 
%     end
% end
% 
% %double line outage, single line outage streams
% E_dbl_single=zeros(length(loseline_dbl), numel(loseline_single)); %E_dbl_single(i,j) gives expected value of log(fj/f0) given double line index i went out 
% for i=1:length(loseline_dbl)
%     if ismember(loseline_dbl(i,:), critical_lines_dbl, 'rows')
%         E_dbl_single(i,:)=NaN;
%         continue
%     end
%     for j=1:numel(loseline_single)
%         if ismember(loseline_single(j), critical_lines_single)
%             E_dbl_single(i,j)=NaN;
%             continue
%         end
%         E_dbl_single(i,j)=KL_compute(sigma_dbl{loseline_dbl(i,1),loseline_dbl(i,2)}, sigma0)-KL_compute(sigma_dbl{loseline_dbl(i,1),loseline_dbl(i,2)}, sigma_single{j}); 
%     end
% end
% 
% %double line outage, double line outage streams
% E_dbl_dbl=zeros(length(loseline_dbl), length(loseline_dbl)); 
% for i=1:length(loseline_dbl)
%     if ismember(loseline_dbl(i,:), critical_lines_dbl, 'rows')
%         E_dbl_dbl(i,:)=NaN;
%         continue
%     end
%     for j=1:length(loseline_dbl)
%         if ismember(loseline_dbl(j,:), critical_lines_dbl, 'rows')
%             E_dbl_dbl(i,j)=NaN;
%             continue
%         end
%         E_dbl_dbl(i,j)=KL_compute(sigma_dbl{loseline_dbl(i,1),loseline_dbl(i,2)}, sigma0)-KL_compute(sigma_dbl{loseline_dbl(i,1),loseline_dbl(i,2)}, sigma_dbl{loseline_dbl(j,1),loseline_dbl(j,2)}); 
%     end
% end
% 
% % save('fourteenbus_dbloutage')
% %%
% %find the largest KL divergences for SINGLE line outages
% index_single_single=cell(length(E_single_single),1);
% index_single_dbl=cell(length(E_single_single),1);
% single_single_mult_pos_drift=zeros(length(E_single_single),1); %num of other streams that also have pos drift due to line i outage
% single_dbl_mult_pos_drift=zeros(length(E_single_single),1); %num of other dbl line streams that also have pos drift due to line i outage
% f_alarm_flag_single=zeros(length(E_single_single),1); %when any other stream's growth rate is larger tha KL divergence of that line
% for i=1:numel(loseline_single)
%     index_single_single{i}=find(E_single_single(i,:)>0);
%     if numel(index_single_single{i})>1
%         single_single_mult_pos_drift(i)=numel(index_single_single{i})-1;
%     end    
%     index_single_dbl{i}=find(E_single_dbl(i,:)>0);
%     if numel(index_single_dbl{i})>0
%         single_dbl_mult_pos_drift(i)=numel(index_single_dbl{i});
%     end  
%     temp=E_single_single(i, 1:end); 
%     temp(i)=[];
%     if KL_single(i)<max([temp E_single_dbl(i,:)])
%         f_alarm_flag_single(i)=1;
%     end
% end
% find(f_alarm_flag_single)
% 
% %%
% %find the lagest KL divergences for DBL line outages
% index_dbl_single=cell(length(E_dbl_single),1);
% index_dbl_dbl=cell(length(E_dbl_single),1);
% dbl_single_mult_pos_drift=zeros(length(E_dbl_single),1); %num of other streams that also have pos drift due to line i outage
% dbl_dbl_mult_pos_drift=zeros(length(E_dbl_single),1); %num of other dbl line streams that also have pos drift due to line i outage
% f_alarm_flag_dbl=zeros(length(E_dbl_single),1); %when any other stream's growth rate is larger tha KL divergence of that line
% 
% for i=1:length(loseline_dbl)
%     index_dbl_single{i}=find(E_dbl_single(i,:)>0);
%     if numel(index_dbl_single{i})>0
%         dbl_single_mult_pos_drift(i)=numel(index_dbl_single{i});
%     end    
%     index_dbl_dbl{i}=find(E_dbl_dbl(i,:)>0);
%     if numel(index_dbl_dbl{i})>1
%         dbl_dbl_mult_pos_drift(i)=numel(index_dbl_dbl{i})-1;
%     end  
%     temp=E_dbl_dbl(i,1:end); 
%     temp(i)=[];
%     if KL_dbl(i)<max([temp E_dbl_single(i,:)])
%         f_alarm_flag_dbl(i)=1;
%     end
% end
% find(f_alarm_flag_dbl)
% 
% 
% count=0; %counts the number of cases where for a double line outage, it's corresponding single line outage streams for the two lines are not highest and second highest
% special_cases=[];
% for i=1:length(loseline_dbl)
%     if ismember(loseline_dbl(i,:), critical_lines_dbl, 'rows')
%         continue
%     end
%     E_dbl_single_sorted=[E_dbl_single(i,:); 1:numel(loseline_single)]';
%     E_dbl_single_sorted=sortrows(E_dbl_single_sorted,-1);
%     temp=isnan(E_dbl_single_sorted);
%     [row col]=find(temp);
%     E_dbl_single_sorted(row,:)=[];
%     if sum(ismember(E_dbl_single_sorted(1:2,2), loseline_dbl(i,:)))~=2
%        count=count+1;
%        if E_dbl_single_sorted(1,1)<14 %may be a good sample path options to produce a case where the two single line outage streams are not the two largest for a dbl outage
%            special_cases=[special_cases; i];
% %            disp(['check dbl outage case ', num2str(i)]);
%        end
%     end
% end
% %percentage of lines where we need to look at dbl line outage streams, doesn't take into account NaN (critical lines)
% loseline_dbl(special_cases,:);
% E_dbl_single(special_cases,:);
% percent=count/length(loseline_dbl) %about 24%
% 
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Useful codes
% 
% % line(1,1:2); %buses of a line
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %given buses, find line
% buses=[4 5]; 
% find((logical(line(:,1)==buses(1))&logical(line(:,2)==buses(2)))|logical((line(:,2)==buses(1))&logical(line(:,1)==buses(2))));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %given lines, find KL div.
% lines_dbl=[linefault_dbl]; 
% index=find((logical(KL_dbl(:,2)==lines_dbl(1))&logical(KL_dbl(:,3)==lines_dbl(2)))|logical((KL_dbl(:,3)==lines_dbl(1))&logical(KL_dbl(:,2)==lines_dbl(2))));
% loseline_dbl(index,:);
% KL_dbl(index,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
