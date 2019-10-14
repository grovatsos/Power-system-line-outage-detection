%%%fourteenbus_DC is the one that runs the 14bus example, generates sample paths of Wn from CuSum algorithm 
%7/16/2014: Automate some of the code so it can be applied to larger
%systems

%7/23/2014: run sample paths for simulations, see second block of this code
%8/5/2014: Note, going from 3 bus system to 14 bus system, the variance on P and Q has been changed,
           %Some lines are better for simulating outages as they affect theta significantly

           
%6/18/2015 creates plots of Detection delay vs mean time to False alarm
%using just mean shift and also CuSum tests

%7/19/2015 Added Y_k variable used to compute markov bound on false alarm
%rates

%CHANGED DELTA THETA SAMPLING from v2



clear all
close all
clc

define_constants;
mpc = case118;
mpopt = mpoption('PF_ALG', 1, 'VERBOSE', 0, 'OUT_ALL', 0);
bus=mpc.bus;
line=mpc.branch;
baseMVA = mpc.baseMVA;
nbus = size(bus,1);
nline = size(line,1);
type = round(bus(:,2));

results = runpf(mpc, mpopt);
theta0=results.bus(:,VA).*pi/180;
theta0(type==3)=[];

% Aggregate multiple lines between the same 2 buses
[u,I,~] = unique(line(:,1:2), 'rows');
ixDupRows = setdiff(1:size(line(:,1:2),1), I);

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

for i = ixDupRows
    from = line(i,1);
    to = line(i,2);
    Xl = line(i,4);
    J0(from,to) = J0(from,to) - 1/Xl;
    J0(to,from) = J0(to,from) - 1/Xl;
end
J0(:,type==3) = []; %remove slack bus
J0(type==3,:) = [];

%Compute single outage KL divergences

%use 0.05 for 14 bus, 0.1 for 3 bus, this does affect KL divergence so in theory, it shouldn't affect ADD
%but, if it hits some limit in power flow and goes to a different operating point the H matrices may change completely 
Psig = 0.05; 
Qsig=0.02;

B0 = inv(J0);
sigma0 = B0*(2*Psig^2*eye(nbus-1))*B0'; %M*Sigma*M' in paper, this is all A1 subsystem

loseline_single = 1:nline;
B_single = cell(numel(loseline_single),1);
sigma= cell(numel(loseline_single),1);
% mu = cell(numel(loseline),1);
% klmean_single = zeros(numel(loseline),1);
KL_single = zeros(numel(loseline_single),2); %KL divergence and associated line loss

%Compute KL divergence for losing a line in A1
for i = 1:numel(loseline_single)
    from = line(loseline_single(i),1);
    to = line(loseline_single(i),2);
    xij = line(loseline_single(i),4);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
    fij(type==3) = []; %remove slack bus
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
    sigma{i} = Ji\(2*Psig^2*eye(nbus-1))/(Ji'); %make it symmetric to avoid numerical rounding issues
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';

    KL_single(i,1) = 1/2 *(trace(sigma0\sigma{i}) - ...
        (length(sigma0)) + log(det(sigma0/sigma{i}))); 
    KL_single(i,2) = loseline_single(i); 
end

   
%Lowest KL divergence dictates the misdetection error. We want KL divergences to be big
KL_single_sorted=sortrows(KL_single,1); %sort KL divergnce from min to max as well as their associated lines
%disp(['min KL divergence for Area 1 single outage = ' num2str(KL_single_sorted(1,1))]);
critical_lines_single=find(isnan(KL_single)); %lines that if lost, creates islanding situations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create sample paths for a single line fault

% seed=4;
% rng(seed);
% linefault=datasample(loseline,1); %choose a random line to fault
linefault_single=74; %(53, 54)
faulttime=Inf;%fault time
if (ismember(linefault_single, critical_lines_single)|| linefault_single<1 || linefault_single > numel(loseline_single))
    error('choose different line to fault')
end


A_single=(100*KL_single_sorted(:,1)/max(KL_single_sorted(:,1)))';
A_single=500*ones(1, numel(loseline_single)); %Threshold


mu0=zeros(1, numel(theta0));
theta0_fault=cell(numel(loseline_single,1));
mu=cell(numel(loseline_single,1));
Y_k=[]; %For markov inequality bound

for i=1:numel(loseline_single)
    if ismember(i, critical_lines_single)
        theta0_fault{i}=NaN;
        continue
    end
    line_temp=line;
    line_temp(loseline_single(i),:)=[]; %remove faulted line
    mpc.branch=line_temp;
    results = runpf(mpc, mpopt);
    theta0_fault_temp=results.bus(:,VA).*pi/180;
    theta0_fault_temp(type==3)=[];
    theta0_fault{i}=theta0_fault_temp;
    mpc.branch=line; %restore original
    mpc.bus=bus;
    mu{i}=theta0_fault{i}-theta0;
end


num_samplepaths=1000;
ADD=Inf;
for pathId=1:num_samplepaths 
    tao_single=1; %time
    f_stop=0;
    f_single=[]; %Specifies which single line outage streams have crossed threshold 
    dyvec = [];
    Wn_Cusum=zeros(1, numel(loseline_single)); %CuSum, one for each line single line outage 
    Wn_Meanshift=zeros(1, numel(loseline_single));
    Wn_Shewhart=zeros(1, numel(loseline_single));
    Wnvec_Cusum = [];
    Wnvec_Meanshift = [];
    Wnvec_Shewhart = [];
    y=theta0; ynew=[]; dy=[];
    %initialize first instant of y. If the simulation starts off with fault, then ynom is faulted y value. 
    if tao_single >= faulttime 
       y=theta0_fault{linefault_single};
    end    
    
    while f_stop==0         
        if tao_single < faulttime %no fault
%             theta_sample= mvnrnd(theta0,sigma0);
%             theta_sample= mvnrnd(zeros(1, numel(theta0)),sigma0);
%             ynew=theta_sample';
%             dy = ynew'-y';
%             dy(type==3)=[];
%             y = ynew;
            dy = mvnrnd(zeros(1, numel(theta0)),sigma0); %Directly generate delta theta from a normal distribution
        else
            
%             theta_sample= mvnrnd(theta0_fault{linefault_single},sigma{linefault_single});
%             ynew=theta_sample';
%             dy = ynew'-y';
%             dy(type==3)=[];
%             y = ynew;
%             dy=mvnrnd(theta0_fault{linefault_single},sigma{linefault_single});
            dy=mvnrnd(zeros(1, numel(theta0)),sigma{linefault_single});
           
        end
        dyvec = [dyvec; dy];

        %Compute the CUSUM statistics for each possible post-change scenario
        for i=1:numel(loseline_single)
            if ismember(loseline_single(i), critical_lines_single)
                continue
            end
%             Wn_Cusum(i) = max(subplus(Wn_Cusum(i) + ...
%             log(mvnpdf(dy, mu0, sigma{i})/mvnpdf(dy, mu0, sigma0))), ...
%             log(mvnpdf(dy,mu{i}',sigma0)/mvnpdf(dy, mu0, sigma0))); 
            
%             Wn_Meanshift(i) = subplus(log(mvnpdf(dy,mu{i}',sigma0)/mvnpdf(dy, mu0, sigma0)));
            

            %SEE COMMENT ABOUT MEAN SHIFT!! 
            Wn_Shewhart(i) = max(log(mvnpdf(dy, mu0, sigma{i})/mvnpdf(dy, mu0, sigma0)), ...
            log(mvnpdf(dy,mu{i}',sigma0)/mvnpdf(dy, mu0, sigma0))); %THIS MEAN SHIFT IS COMPUTED USING OLD METHOD (through power flow), But we sample delta theta directly from the Norm distribution with 0 mean
        end
%         Wnvec_Cusum = [Wnvec_Cusum; Wn_Cusum];
%         Wnvec_Meanshift = [Wnvec_Meanshift; Wn_Meanshift];
        Wnvec_Shewhart = [Wnvec_Shewhart; Wn_Shewhart];
       
%         if sum(subplus(Wn_Cusum-A_single))>0
%       if sum(subplus(Wn_Meanshift-A_single))>0
      
      Y_k=[Y_k max(Wn_Shewhart)];

      if sum(subplus(Wn_Shewhart-A_single))>0
         
            [max_temp max_single]=max(Wn_Shewhart-A_single);
            Sorted_single=sort(subplus(Wn_Shewhart-A_single), 'descend');
%               find(temp(1)==(Wn_single-A_single))
            f_single=find(subplus(Wn_Shewhart-A_single)); 
            disp(['single line outage streams ', num2str(f_single), ' crossed the treshold ', num2str(A_single(f_single))]);
            
%             if tao_single<faulttime || max_single~=linefault_single
%                error('misdetection!')
%             end
            
            if pathId==1 
                disp(['tao_single = ' num2str(tao_single)]);
%                 ADD = tao_single-faulttime+1;
                ADD = tao_single %For false alarm detection rate computation since fault time is Inf
            else
                disp(['tao_single = ' num2str(tao_single)]);
%                 ADD = ((pathId-1)*ADD + (tao_single-faulttime+1))/pathId; %running average to compute ADD over many sample paths
                ADD = ((pathId-1)*ADD + tao_single)/pathId %For false alarm detection rate computation since fault time is Inf
            end
            pathId;
            f_stop==1;
            break; %needed so tao_single doesn't increase  
        end 
        %Update the stopping variable till we actually stop
        if mod(tao_single, 200)==0
            disp(['tao_single = ', num2str(tao_single)]);
            Y_k;
        end
        
        tao_single=tao_single+1;
       
        
        if tao_single>5000
            mean(Y_k)
            error('timed out')
        end    
        
    end
pathId    
end
disp(['single line ' num2str(linefault_single), ' fault ADD=' num2str(ADD)]);

figure('position', [20 80 1300 800]);
axes('position', [0.14 0.14 0.8 0.8]);
hold on; box on; grid on;
title(['Single Line Outage ', num2str(linefault_single)])
p1=plot(1:tao_single, Wnvec_Shewhart(1:tao_single, linefault_single)','--b', 'LineWidth', 2, 'Marker', '*', 'Markersize', 6);
p2=plot(1:tao_single, Wnvec_Shewhart(1:tao_single,5)','--g', 'LineWidth', 2,'Marker', 'o','Markersize', 6);
p3=plot(1:tao_single, Wnvec_Shewhart(1:tao_single,2)','--r', 'LineWidth', 2,'Marker', 's','Markersize', 6);
p4=plot(1:tao_single, Wnvec_Shewhart(1:tao_single,6)','--c', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% p5=plot(1:tao_single, Wnvec_dbl(1:tao_single,max_dbl)','--m', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% l1=legend([p1, p2, p3, p4, p5], '$W_{linefault single}$','$W_{5}$','$W_{2}$', '$W_{6}$', '$W_{max dbl}$');
% p1=plot(1:tao_single, Wnvec_Meanshift(1:tao_single, linefault_single)','-b', 'LineWidth', 2, 'Marker', '*', 'Markersize', 6);

set(gca,'fontsize',36,'fontname','times new roman'); ...
   xl=xlabel('$k$','fontsize',38,'fontname','times new roman');
   set(xl,'Interpreter','latex');
   yl=ylabel('$W[k]$','fontsize',38,'fontname','times new roman');
   set(yl,'Interpreter','latex');
%    set(l1,'Interpreter','latex');
%    set(l1,'FontSize',37, 'Interpreter','latex')
%    set(gca,'XTick',1.0487:0.0001:1.0492)
%    set(gca,'YTick',4952270:40:4952290)
%    axis([1.04901 1.0504 0.98 1.07])

figure
grid on
plot(1:tao_single, Wnvec_Shewhart(1:tao_single,:)');






