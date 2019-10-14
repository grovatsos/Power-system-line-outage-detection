% This is the code for simulating the shewhart test, cusum, meanshift, covariance shift. It samples from a
% prechange and post change distribution. 
% Simulates the avergage ADD
%8/30/2015

%9/7/2015: Modified this code to only simulate ADD for meanshift. Instead of
%sampling from covariance shift distribution for importance sampling, we
%sample from a meanshift distribution. USES importance sampling (see notes
%for details)


clear all
close all
% clc

define_constants;
mpc = case14;
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
Psig = 0.5; 
Qsig=0.01;

B0 = inv(J0);
sigma0 = B0*(2*Psig^2*eye(nbus-1))*B0'; %M*Sigma*M' in paper, this is all A1 subsystem

loseline = 1:nline;
B_single = cell(numel(loseline),1);
sigma= cell(numel(loseline),1);
% mu = cell(numel(loseline),1);
% klmean_single = zeros(numel(loseline),1);
KL_single = zeros(numel(loseline),2); %KL divergence and associated line loss

%Compute KL divergence for losing a line in A1
for i = 1:numel(loseline)
    from = line(loseline(i),1);
    to = line(loseline(i),2);
    xij = line(loseline(i),4);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
    fij(type==3) = []; %remove slack bus
    dJ = -1/xij*(fij*fij');
    Ji = J0 + dJ; %H matrix in paper
    
    if rank(Ji) < size((Ji),1)
        disp(['critical single line ' num2str(loseline(i))]);
        disp(['line ' num2str(loseline(i)), ' has buses ' num2str(from), ' and ' num2str(to), sprintf('\n')]);
        B_single{i}=NaN;
        sigma_single{i}=NaN;
        KL_single(i,1) = NaN; 
        KL_single(i,2) = loseline(i); 

        continue
    end
    
    B_single{i} = inv(Ji); %Compute each M corresponding to each line outage in A1;
    sigma{i} = Ji\(2*Psig^2*eye(nbus-1))/(Ji'); %make it symmetric to avoid numerical rounding issues
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';

    KL_single(i,1) = 1/2 *(trace(sigma0\sigma{i}) - ...
        (length(sigma0)) + log(det(sigma0/sigma{i}))); 
    KL_single(i,2) = loseline(i); 
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
linefault=7; %used for importance sampling
faulttime=0;%fault time


% A=100*ones(1, numel(loseline)); %Threshold

%%
%gamma=10000;
%A=nline*gamma*ones(1,numel(loseline));

if (ismember(linefault, critical_lines_single)|| linefault<1 || linefault > numel(loseline))
    error('choose different line to fault')
end


% A_single=(100*KL_single_sorted(:,1)/max(KL_single_sorted(:,1)))';
% A_single=175*ones(1, numel(loseline)); %Threshold

% %"Optimal" Threshold from George's derivation
% % FA_rate=120;
% zeta_theta=0.0274; %obtained from saved .mat files and George_threshold_v2 code
% gamma=5e2;
% alpha=log(gamma)+log(numel(loseline));
% A_single=alpha-log(1./((KL_single(:,1))'.*zeta_theta.^2));


mu0=zeros(1, numel(theta0));
theta0_fault=cell(numel(loseline,1));
mu_fault=cell(numel(loseline,1));

for i=1:numel(loseline)
    if ismember(i, critical_lines_single)
        theta0_fault{i}=NaN;
        continue
    end
    line_temp=line;
    line_temp(loseline(i),:)=[]; %remove faulted line
    mpc.branch=line_temp;
    results = runpf(mpc, mpopt);
    theta0_fault_temp=results.bus(:,VA).*pi/180;
    theta0_fault_temp(type==3)=[];
    theta0_fault{i}=theta0_fault_temp;
    mpc.branch=line; %restore original
    mpc.bus=bus;
    mu_fault{i}=(theta0_fault{i}-theta0)';
end

%%
num_samplepaths=1;
num_samples=5000;

threshold=[1 2 3 4 5 6 7]; 


Prob_vec_mean=zeros(numel(threshold), num_samplepaths);
Prob_vec_cov=zeros(numel(threshold), num_samplepaths);
ADD_vec = zeros(numel(threshold), num_samplepaths);

for pathId=1:num_samplepaths 
%     seed=pathId;
%     rng(seed);
    Y_k_mean=[]; 
%     Y_k_cov=[];
    tao=1; %time
    f_stop=0;
%     f_single=[]; %Specifies which single line outage streams have crossed threshold 
    dyvec_mean = [];
    dyvec_cov = [];
    Wn_Cusum=zeros(1, numel(loseline)); %CuSum, one for each line single line outage 
    Wn_Meanshift=zeros(1, numel(loseline));
    Wn_Shewhart=zeros(1, numel(loseline));
    Wn_Covshift=zeros(1, numel(loseline));
    Wnvec_Cusum = [];
    Wnvec_Meanshift = [];
    Wnvec_Shewhart = [];
    Wnvec_Covshift = [];
    y=theta0; ynew=[]; dy_mean=[]; dy_cov=[];
    %initialize first instant of y. If the simulation starts off with fault, then ynom is faulted y value. 
%     if tao >= faulttime 
%        y=theta0_fault{linefault};
%     end    
    
    while f_stop==0         
        if tao < faulttime %no fault
%             theta_sample= mvnrnd(zeros(1, numel(theta0)),sigma0);
%             ynew=theta_sample';
%             dy = ynew'-y';
%             dy(type==3)=[];
%             y = ynew;
            dy = mvnrnd(zeros(1, numel(theta0)),sigma0);
        elseif tao ==faulttime %at the fault instance
%             theta_sample= mvnrnd(mu{linefault_single},sigma0);
%             ynew=theta_sample';
%             dy = ynew'-y';
%             dy(type==3)=[];
%             y = ynew;
            dy = mvnrnd(mu_fault{linefault},sigma0);
        else
%             theta_sample= mvnrnd(theta0_fault{linefault_single},sigma{linefault_single}); %check mean of distribution. 
%             ynew=theta_sample';
%             dy = ynew'-y';
%             dy(type==3)=[];
%             y = ynew;
%             dy_cov=mvnrnd(zeros(1, numel(theta0)),sigma{linefault});            
            
            %%%%%%%%%%%%%%%%%%%%%
            dy_mean=mvnrnd(mu_fault{linefault},sigma0); %Sample from meanshift distribution for importance sampling
            %%%%%%%%%%%%%%%%%%%%%%%
           
        end
        dyvec_mean = [dyvec_mean; dy_mean];
%         dyvec_cov = [dyvec_cov; dy_cov];

        %Compute the CUSUM statistics for each possible post-change scenario
        for i=1:numel(loseline)
%             if ismember(loseline_single(i), critical_lines_single)
%                 continue
%             end
%             Wn_Cusum(i) = max(subplus(Wn_Cusum(i) + ...
%             log(mvnpdf(dy, mu0, sigma{i})/mvnpdf(dy, mu0, sigma0))), ...
%             log(mvnpdf(dy,mu_fault{i},sigma0)/mvnpdf(dy, mu0, sigma0))); 
            
            Wn_Meanshift(i) = subplus(log(mvnpdf(dy_mean,mu_fault{i},sigma0)/mvnpdf(dy_mean, mu0, sigma0)));
%             Wn_Covshift(i) = subplus(log(mvnpdf(dy_cov, mu0, sigma{i})/mvnpdf(dy_cov, mu0, sigma0)));
            
%             Wn_Shewhart(i) = max(log(mvnpdf(dy, mu0, sigma{i})/mvnpdf(dy, mu0, sigma0)),...
%             log(mvnpdf(dy,mu_fault{i},sigma0)/mvnpdf(dy, mu0, sigma0)));
        end
%         Wnvec_Cusum = [Wnvec_Cusum; Wn_Cusum];
        Wnvec_Meanshift = [Wnvec_Meanshift; Wn_Meanshift];
%         Wnvec_Shewhart = [Wnvec_Shewhart; Wn_Shewhart];
%         Wnvec_Covshift = [Wnvec_Covshift; Wn_Covshift];
       
       
      
%       Y_k=[Y_k max(Wn_Shewhart)];
        Y_k_mean=[Y_k_mean Wn_Meanshift(linefault)];
%         Y_k_cov=[Y_k_cov max(Wn_Covshift)];
 

%         if sum(subplus(Wn_Cusum-A))>0
% %       if sum(subplus(Wn_Meanshift-A))>0
% %       if sum(subplus(Wn_Shewhart-A))>0
% %         if sum(subplus(Wn_Covshift-A))>0
%          
%             [max_temp max_single]=max(Wn_Cusum-A);
%             Sorted_single=sort(subplus(Wn_Cusum-A), 'descend');
% %               find(temp(1)==(Wn_single-A_single))
%             f_single=find(subplus(Wn_Cusum-A)); 
%             
%             if (~ismember(linefault, f_single))
%                 error('misdetection occured during importance sampling');
%             end
%             
%             disp(['Single line outage streams ', num2str(f_single), ' crossed the treshold ', num2str(A(f_single))]);
%             
% %             if tao_single<faulttime || max_single~=linefault_single
% %                error('misdetection!')
% %             end
%             
%             if pathId==1 
%                 disp(['tao = ' num2str(tao)]);
%                 FA = tao;
%             else
%                 disp(['tao = ' num2str(tao)]);
% %                 ADD = ((pathId-1)*ADD + (tao_single-faulttime+1))/pathId; %running average to compute ADD over many sample paths
%                 FA = ((pathId-1)*FA + tao)/pathId %For false alarm detection rate computation since fault time is Inf
%             end
%             pathId;
%             f_stop=1;
%             break; %needed so tao_single doesn't increase  
%         end 
        %Update the stopping variable till we actually stop
        if mod(tao, 2000)==0
            disp(['tao = ', num2str(tao)]);
%             mean(Y_k)
        end
        
        tao=tao+1;
        

        
        if tao>num_samples
            f_stop=1;
%             mean(Y_k)
%             error('timed out')
        end    
        
    end
    if mod(pathId, 3000)==0
            percent=pathId/num_samples*100;
            disp(['percent done = ', num2str(percent)]);
%             mean(Y_k)
    end
    
    
    
for i=1:numel(threshold) 
count_mean=0;
temp=0;
    for j=1:numel(Y_k_mean)
        if Y_k_mean(j)>=threshold(i)
           temp=temp+mvnpdf(dyvec_mean(j,:), mu0, sigma{linefault})/mvnpdf(dyvec_mean(j,:), mu_fault{linefault},sigma0); %importance sampling
           count_mean=count_mean+1;
        end    
        
    end
    Prob_vec_mean(i, pathId)=count_mean/numel(Y_k_mean);
    Prob_vec_cov(i, pathId)=temp/numel(Y_k_mean);   
    ADD_vec(i, pathId)=Prob_vec_mean(i, pathId)+(1+1./Prob_vec_cov(i, pathId))*(1-Prob_vec_mean(i, pathId));
    i
end
end  
threshold
ADD_vec %each row corresponds to a threshold 

