%%%fourteenbus_DC is the one that runs the 14bus example, generates sample paths of Wn from CuSum algorithm 
%7/16/2014: Automate some of the code so it can be applied to larger
%systems
%7/23/2014: run sample paths for simulations, see second block of this code
%8/5/2014: Note, going from 3 bus system to 14 bus system, the variance on P and Q has been changed,
           %Some lines are better for simulating outages as they affect theta significantly
%10/28/2014: This code checks ADD by setting individual thresholds for each line outage scaled according to their KL divergence
%11/18/2014: Adapted from 3bus to 14bus system. 
%3/17/2015: Adapted to include double line outages
%4/28/2015: Uses distributed slack bus
%uses Matpower to generate sample plots
% This code allows for PMUs at a subset of buses in addition to computing sample paths

clear all
% close all
% clc

%% Initialize case
define_constants;
mpc = case14;
mpopt = mpoption('PF_ALG', 1, 'VERBOSE', 0, 'OUT_ALL', 0);
bus=mpc.bus;
line=mpc.branch;
baseMVA = mpc.baseMVA;
nbus = size(bus,1);
nline = size(line,1);
type = round(bus(:,2));
genbus=(type~=1);
genbus_noslack=(type==2);
loadbus=(type==1);
slackbus=(type==3);


%create participation factor matrix, scale according to Pmax
alpha=ones(numel(find(genbus)),numel(find(loadbus)));
[l,~]=size(alpha);
for i=1:l
    alpha(i,:)=mpc.gen(i,9)/sum(mpc.gen(:,9));
end

alpha(1,:)=[]; %remove slack bus

results = runpf(mpc, mpopt);
theta0=results.bus(:,VA).*pi/180;

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
J0_orig=J0;
J0(:,type==3) = []; %remove slack bus
J0(type==3,:) = [];

% create the new jacobian matrix based on new bus indices [Pg Pd]= J*[theta_gen theta_load]
J0_dsb=J0_orig*0;

gen_indices=find(genbus_noslack);
load_indices=find(loadbus);
slack_index=find(slackbus);
for i=1:length(J0_dsb)
    for j=1:length(J0_dsb)
        J0_ij=J0_orig(i,j)   
        i_new=dsb_bus(i, gen_indices, load_indices, slack_index);
        j_new=dsb_bus(j, gen_indices, load_indices, slack_index);
        J0_dsb(i_new,j_new)=J0_ij;    
    end
end

J0_dsb(:,1) = []; %remove slack bus
J0_dsb(1,:) = [];

%Compute single outage KL divergences

%use 0.05 for 14 bus, 0.1 for 3 bus, this does affect KL divergence so in theory, it shouldn't affect ADD
%but, if it hits some limit in power flow and goes to a different operating point the H matrices may change completely 
Psig = 0.03; 
Qsig=0.02;

B0_dsb = inv(J0_dsb); %M0 in paper
M0_1=B0_dsb(:,1:numel(gen_indices));
M0_2=B0_dsb(:,(numel(gen_indices)+1):end);
M0_hat=M0_1*alpha+M0_2;

%rearrage the rows of M0_hat so that the buses for theta corresponds to original bus indices (simplifies the code later on)
M0_hat=rearrage_M_hat(M0_hat,gen_indices,load_indices,slack_index);

remove_PMU_bus=[1 2 3 4 5 7]; %buses that do not have PMUs, note that here, slack has already been removed so bus1 is the first generator bus that's not slack bus must remove enough buses s.t. dim of theta=dim of Pd
% remove_PMU_bus=[];
if (nbus-1-numel(remove_PMU_bus))>numel(load_indices)
    error('remove more buses with PMU measurements');
end

% Compute the KL divergences with certain bus' PMUs removed
PMU_bus=setdiff([1:(nbus-1)], remove_PMU_bus); %bus that has the PMU
M0_hat(remove_PMU_bus,:)=[];
sigma0 = M0_hat*(2*Psig^2*eye(numel(load_indices)))*M0_hat'; %Note that Covariance matrix is smaller due to only load buses being stimulated

loseline_single = 1:nline;
B_single = cell(numel(loseline_single),1);
sigma_single = cell(numel(loseline_single),1);
% mu = cell(numel(loseline),1);
% klmean_single = zeros(numel(loseline),1);
KL_single = zeros(numel(loseline_single),2); %KL divergence and associated line loss

%Compute KL divergence for losing a line in A1
for i = 1:numel(loseline_single)
    from = line(loseline_single(i),1);
    to = line(loseline_single(i),2);
    xij = line(loseline_single(i),4);
    fij = zeros(nbus,1);
    from_new=dsb_bus(from, gen_indices, load_indices, slack_index);
    to_new=dsb_bus(to, gen_indices, load_indices, slack_index);
    fij(from_new) = 1;
    fij(to_new) = -1;
    fij(1) = []; %remove slack bus
    dJ_dsb = -1/xij*(fij*fij');
    Ji_dsb = J0_dsb + dJ_dsb; %H matrix in paper
    
    if rank(Ji_dsb) < size((Ji_dsb),1)
        disp(['critical single line ' num2str(loseline_single(i))]);
        disp(['line ' num2str(loseline_single(i)), ' has buses ' num2str(from), ' and ' num2str(to), sprintf('\n')]);
        B_single{i}=NaN;
        sigma_single{i}=NaN;
        KL_single(i,1) = NaN; 
        KL_single(i,2) = loseline_single(i); 

        continue
    end
    
    Bi_dsb=inv(Ji_dsb); %Compute each M corresponding to each line outage in A1;
    Mi_1=Bi_dsb(:,1:numel(gen_indices));
    Mi_2=Bi_dsb(:,(numel(gen_indices)+1):end);
    Mi_hat=Mi_1*alpha+Mi_2;
    
    Mi_hat=rearrage_M_hat(Mi_hat,gen_indices,load_indices,slack_index);
    
    
    Mi_hat(remove_PMU_bus,:)=[]; %take out buses without PMU measurements

    sigma_single{i} = Mi_hat*(2*Psig^2*eye(numel(load_indices)))*Mi_hat'; %make it symmetric to avoid numerical rounding issues
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



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create sample paths for a single line fault

seed=4;
% rng(seed);
% linefault=datasample(loseline,1); %choose a random line to fault
linefault_single=10;
faulttime=1;%fault time
if (ismember(linefault_single, critical_lines_single)|| linefault_single<1 || linefault_single > numel(loseline_single))
    error('choose different line to fault')
end
num_samplepaths=1;
y=theta0;
ADD=Inf;
%run sample paths
%scale thresholds according to KL divergences
% A_single=(100*KL_single_limitedPMUs(:,1)/max(KL_single_limitedPMUs(:,1)))';
A_single=100*ones(1, numel(loseline_single)); %Threshold


mu0=zeros(1, numel(theta0));
theta0_fault=cell(numel(loseline_single,1));
for i=1:numel(loseline_single)
    if ismember(i, critical_lines_single)
        theta0_fault{i}=NaN;
        continue
    end
    line_temp=line;
    line_temp(loseline_single(i),:)=[]; %remove faulted line
    mpc.branch=line_temp;
    results = runpf(mpc, mpopt);
    theta0_fault{i}=results.bus(:,VA).*pi/180;
    mpc.branch=line; %restore original
    mpc.bus=bus;
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
        Pd=bus(type==1,3); %Only load buses generate powerNegative signs are to match earlier version of code. Note that Pload and Qload may become negative meaning bus is generating power. 
        Pg=bus(type==2,3);
        Q=bus(type~=3,4);
        delta_Pd=randn(numel(Pd),1)*Psig*baseMVA;
        Pd_sample=Pd+delta_Pd;
        Pg_sample=Pg+alpha*delta_Pd; %governor power flow
        Q_sample=Q+randn(numel(Q),1)*Qsig*baseMVA;
         
        if tao_single < faulttime %no fault
            bus_sample=bus;
            bus_sample(type==1,3)=Pd_sample;
            bus_sample(type==2,3)=Pg_sample;
            bus_sample(type~=3,4)=Q_sample;
            mpc.bus=bus_sample;
            results = runpf(mpc, mpopt);
            theta_sample=results.bus(:,VA).*pi/180;
            mpc.branch=line; %restore original
            mpc.bus=bus;
            ynew=theta_sample;
            dy = ynew'-y';
            dy(type==3)=[];
            y = ynew;
        else
            bus_sample=bus;
            bus_sample(type==1,3)=Pd_sample;
            bus_sample(type==2,3)=Pg_sample;
            bus_sample(type~=3,4)=Q_sample;
            line_temp=line;
            line_temp(linefault_single,:)=[];
            mpc.bus=bus_sample;
            mpc.branch=line_temp;
            results = runpf(mpc, mpopt);
            theta_sample=results.bus(:,VA).*pi/180;
            mpc.branch=line; %restore original
            mpc.bus=bus;
            ynew=theta_sample;
            dy = ynew'-y';
            dy(type==3)=[];
            y = ynew;
           
        end
        dyvec = [dyvec; dy];

        %Compute the CUSUM statistics for each possible post-change scenario
        for i=1:numel(loseline_single)
            if ismember(loseline_single(i), critical_lines_single)
                continue
            end
            Wn_single(i) = subplus(Wn_single(i) + log(mvnpdf(dy(PMU_bus), mu0(PMU_bus), sigma_single{i})/mvnpdf(dy(PMU_bus), mu0(PMU_bus), sigma0)));
        end
        Wnvec_single = [Wnvec_single; Wn_single];
        

            
        
        %Check if any one of them is above A
        if sum(subplus(Wn_single-A_single))>0
            [temp_single, max_single]=max(Wn_single-A_single);
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
            f_stop=1;
            break; %needed so tao_single doesn't increase  
        end 
        %Update the stopping variable till we actually stop
        tao_single=tao_single+1;
        
        if tao_single>1000
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
p1=plot(1:tao_single, Wnvec_single(1:tao_single, linefault_single)','--b', 'LineWidth', 2, 'Marker', '*', 'Markersize', 6);
p2=plot(1:tao_single, Wnvec_single(1:tao_single,4)','--g', 'LineWidth', 2,'Marker', 'o','Markersize', 6);
p3=plot(1:tao_single, Wnvec_single(1:tao_single,5)','--r', 'LineWidth', 2,'Marker', 's','Markersize', 6);
p4=plot(1:tao_single, Wnvec_single(1:tao_single,6)','--c', 'LineWidth', 2,'Marker', 's','Markersize', 6);
l1=legend([p1, p2, p3, p4], '$W_{linefault single}$','$W_{4}$','$W_{5}$', '$W_{6}$');

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
% pause

