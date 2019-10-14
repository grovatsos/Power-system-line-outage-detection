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

%7/27/2015: Use imporatance sampling technique, see 7/26/2015 notes in
%notebook

%8/20/2015: Compute the variance of the samples as well so we know how good
%our estimate is 

%8/23/2015: Modified the code so that it generates the plots of Shewhart
%test and Cusum tests. Compute FA. Also, the meanshift is computed
%from nonlinear power flow but the samples of theta are taken directly from
%normal distribution. 

%8/24/2015 Use importance sampling technique. It only works for the
%Shewhard test since it's a one-shot detection. Not sure how to work it in
%for the CuSum test but the FA for CuSums can be found by computing the
%E_f0(log(f_i/f_0)) and the variance of log(f_i/f_0) under f_0 distribution









clear all
close all
% clc

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
Psig = 0.5; 
Qsig=0.01;

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
linefault_single=7; %used for importance sampling
faulttime=1;%fault time
if (ismember(linefault_single, critical_lines_single)|| linefault_single<1 || linefault_single > numel(loseline_single))
    error('choose different line to fault')
end


% A_single=(100*KL_single_sorted(:,1)/max(KL_single_sorted(:,1)))';
% A_single=175*ones(1, numel(loseline_single)); %Threshold

% %"Optimal" Threshold from George's derivation
% % FA_rate=120;
% zeta_theta=0.0274; %obtained from saved .mat files and George_threshold_v2 code
% gamma=5e2;
% alpha=log(gamma)+log(numel(loseline_single));
% A_single=alpha-log(1./((KL_single(:,1))'.*zeta_theta.^2));


mu0=zeros(1, numel(theta0));
theta0_fault=cell(numel(loseline_single,1));
mu=cell(numel(loseline_single,1));
% Y_k=[]; 


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


num_samplepaths=1;
num_samples=300000;
% FA=Inf;
threshold=[1 2 3 4 5]; %
FA_prob_vec=zeros(numel(threshold), num_samplepaths);
FA_vec=zeros(numel(threshold), num_samplepaths);

for pathId=1:num_samplepaths 
%     seed=pathId;
%     rng(seed);
    Y_k=[];
    tao_single=1; %time
    f_stop=0;
    f_single=[]; %Specifies which single line outage streams have crossed threshold 
    dyvec = [];
    Wn_Cusum=zeros(1, numel(loseline_single)); %CuSum, one for each line single line outage 
    Wn_Meanshift=zeros(1, numel(loseline_single));
    Wn_Shewhart=zeros(1, numel(loseline_single));
    Wn_Covshift=zeros(1, numel(loseline_single));
    Wnvec_Cusum = [];
    Wnvec_Meanshift = [];
    Wnvec_Shewhart = [];
    Wnvec_Covshift = [];
    y=theta0; ynew=[]; dy=[];
    %initialize first instant of y. If the simulation starts off with fault, then ynom is faulted y value. 
    if tao_single >= faulttime 
       y=theta0_fault{linefault_single};
    end    
    
    while f_stop==0         
        if tao_single < faulttime %no fault
%             theta_sample= mvnrnd(zeros(1, numel(theta0)),sigma0);
%             ynew=theta_sample';
%             dy = ynew'-y';
%             dy(type==3)=[];
%             y = ynew;
            dy = mvnrnd(zeros(1, numel(theta0)),sigma0);
        elseif tao_single ==faulttime %at the fault instance
%             theta_sample= mvnrnd(mu{linefault_single},sigma0);
%             ynew=theta_sample';
%             dy = ynew'-y';
%             dy(type==3)=[];
%             y = ynew;
            dy = mvnrnd(mu{linefault_single},sigma0);
        else
%             theta_sample= mvnrnd(theta0_fault{linefault_single},sigma{linefault_single}); %check mean of distribution. 
%             ynew=theta_sample';
%             dy = ynew'-y';
%             dy(type==3)=[];
%             y = ynew;
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
%             Wn_Covshift(i) = subplus(log(mvnpdf(dy, mu0, sigma{i})/mvnpdf(dy, mu0, sigma0)));
            
            Wn_Shewhart(i) = max(log(mvnpdf(dy, mu0, sigma{i})/mvnpdf(dy, mu0, sigma0)),...
            log(mvnpdf(dy,mu{i}',sigma0)/mvnpdf(dy, mu0, sigma0)));
        end
%         Wnvec_Cusum = [Wnvec_Cusum; Wn_Cusum];
%         Wnvec_Meanshift = [Wnvec_Meanshift; Wn_Meanshift];
        Wnvec_Shewhart = [Wnvec_Shewhart; Wn_Shewhart];
%         Wnvec_Covshift = [Wnvec_Covshift; Wn_Covshift];
       
       
      
      Y_k=[Y_k max(Wn_Shewhart)];
%          Y_k=[Y_k max(Wn_Covshift)];

%         if sum(subplus(Wn_Cusum-A_single))>0
% %       if sum(subplus(Wn_Meanshift-A_single))>0
% %       if sum(subplus(Wn_Shewhart-A_single))>0
% %         if sum(subplus(Wn_Covshift-A_single))>0
%          
%             [max_temp max_single]=max(Wn_Cusum-A_single);
%             Sorted_single=sort(subplus(Wn_Cusum-A_single), 'descend');
% %               find(temp(1)==(Wn_single-A_single))
%             f_single=find(subplus(Wn_Cusum-A_single)); 
%             disp(['Single line outage streams ', num2str(f_single), ' crossed the treshold ', num2str(A_single(f_single))]);
%             
% %             if tao_single<faulttime || max_single~=linefault_single
% %                error('misdetection!')
% %             end
%             
%             if pathId==1 
%                 disp(['tao_single = ' num2str(tao_single)]);
%                 FA = (tao_single-faulttime+1)*min((mvnpdf(dy, mu0, sigma{linefault_single})/mvnpdf(dy, mu0, sigma0)),...
%                     (mvnpdf(dy, mu{linefault_single}', sigma0)/mvnpdf(dy, mu0, sigma0)));
%             else
%                 disp(['tao_single = ' num2str(tao_single)]);
%                 FA_temp=(tao_single-faulttime+1)*min((mvnpdf(dy, mu0, sigma{linefault_single})/mvnpdf(dy, mu0, sigma0)),...
%                         (mvnpdf(dy, mu{linefault_single}', sigma0)/mvnpdf(dy, mu0, sigma0)));
%                 FA = (((pathId-1)*FA + FA_temp)/pathId); %running average to compute ADD over many sample paths
% %                 FA = ((pathId-1)*FA + tao_single)/pathId %For false alarm detection rate computation since fault time is Inf
%             end
%             pathId;
%             f_stop=1;
%             break; %needed so tao_single doesn't increase  
%         end 
        %Update the stopping variable till we actually stop
        if mod(tao_single, 2000)==0
            disp(['tao_single = ', num2str(tao_single)]);
%             mean(Y_k)
        end
        
        tao_single=tao_single+1;
        

        
        if tao_single>num_samples
            f_stop=1;
%             mean(Y_k)
%             error('timed out')
        end    
        
    end
%     if mod(pathId, 20)==0
%             disp(['current path = ', num2str(pathId)]);
% %             mean(Y_k)
%     end

% threshold=[75 100]; %
% FA_prob_vec=zeros(numel(threshold), num_samplepaths);
% FA_vec=zeros(numel(threshold), num_samplepaths);

for i=1:numel(threshold) 
temp=0;
count=0;
    for j=1:numel(Y_k)
        if Y_k(j)>=threshold(i)
            temp=temp+mvnpdf(dyvec(j,:), mu0, sigma0)/mvnpdf(dyvec(j,:), mu0, sigma{linefault_single});
            count=count+1;
        end
%         E_l(i)=temp/numel(Y_k);
    end
    FA_prob_vec(i, pathId)=temp/numel(Y_k);
    FA_vec(i, pathId)=1./FA_prob_vec(i, pathId);
end
end

FA_vec




% 








%%


% disp(['single line ' num2str(linefault_single), ' fault ADD=' num2str(ADD)]);
% 
% figure('position', [20 80 1300 800]);
% axes('position', [0.14 0.14 0.8 0.8]);
% hold on; box on; grid on;
% title(['Single Line Outage ', num2str(linefault_single)])
% p1=plot(1:tao_single, Wnvec_Cusum(1:tao_single, linefault_single)','--b', 'LineWidth', 2, 'Marker', '*', 'Markersize', 6);
% p2=plot(1:tao_single, Wnvec_Cusum(1:tao_single,5)','--g', 'LineWidth', 2,'Marker', 'o','Markersize', 6);
% p3=plot(1:tao_single, Wnvec_Cusum(1:tao_single,2)','--r', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% p4=plot(1:tao_single, Wnvec_Cusum(1:tao_single,6)','--c', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% % p5=plot(1:tao_single, Wnvec_dbl(1:tao_single,max_dbl)','--m', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% % l1=legend([p1, p2, p3, p4, p5], '$W_{linefault single}$','$W_{5}$','$W_{2}$', '$W_{6}$', '$W_{max dbl}$');
% % p1=plot(1:tao_single, Wnvec_Meanshift(1:tao_single, linefault_single)','-b', 'LineWidth', 2, 'Marker', '*', 'Markersize', 6);
% 
% set(gca,'fontsize',36,'fontname','times new roman'); ...
%    xl=xlabel('$k$','fontsize',38,'fontname','times new roman');
%    set(xl,'Interpreter','latex');
%    yl=ylabel('$W[k]$','fontsize',38,'fontname','times new roman');
%    set(yl,'Interpreter','latex');
% %    set(l1,'Interpreter','latex');
% %    set(l1,'FontSize',37, 'Interpreter','latex')
% %    set(gca,'XTick',1.0487:0.0001:1.0492)
% %    set(gca,'YTick',4952270:40:4952290)
% %    axis([1.04901 1.0504 0.98 1.07])
% 
% figure
% grid on
% plot(1:tao_single, Wnvec_Shewhart(1:tao_single,:)');
% 
% 
% 
% 
% % %Compute KL divergences and also some False Isolation parameters
% % KL1=KL_compute(sigma{1}, sigma0); 
% % KL2=KL_compute(sigma{2}, sigma0); 
% % KL3=KL_compute(sigma{3}, sigma0);
% % KL_fault=KL_compute(sigma{linefault}, sigma0);
% % E120=1/2*(trace(sigma0\sigma{1})-trace(sigma{2}\sigma{1})+ log(det(sigma0/sigma{2}))); 
% % E130=1/2*(trace(sigma0\sigma{1})-trace(sigma{3}\sigma{1})+ log(det(sigma0/sigma{3})));
% % E210=KL_compute(sigma2, sigma0)-KL_compute(sigma2, sigma1); 
% % % E310=1/2*(trace(sigma0\sigma3)-trace(sigma1\sigma3)+ log(det(sigma0/sigma1))); 
% % % E320=1/2*(trace(sigma0\sigma3)-trace(sigma2\sigma3)+ log(det(sigma0/sigma2))); 
% 
