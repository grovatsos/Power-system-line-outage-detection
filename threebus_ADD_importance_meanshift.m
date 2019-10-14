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
clc


[bus, line] = buslinedata(3);
nbus = size(bus,1);
nline = size(line,1);
type = round(bus(:,10));

[soln, ybus, J] = get_pf(bus,line);
theta0 = soln(1:end,3).*pi/180; %in radians
theta0(type==1)=[];

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


J0_orig=J0;
J0(:,type==1) = [];   J0(type==1,:) = []; %remove slack bus

Psig = 0.5;
B0=inv(J0);
loseline = 1:nline;


mu0 = zeros(1, nbus-1); %don't include slack
sigma0 = B0*(2*Psig^2*eye(nbus-1))*B0'; %M*Sigma*M' in paper

B = cell(numel(loseline),1);
sigma = cell(numel(loseline),1);
mu_fault = cell(numel(loseline),1);
theta0_fault=cell(numel(loseline,1));
klcov_lines = zeros(numel(loseline),3); %KL divergence and associated line loss


for i = 1:numel(loseline)
    
    line2=line;
    line2(i,:)=[]; %remove faulted line
    [soln_fault, ybus_fault, J_fault] = get_pf(bus,line2);
    theta0_fault_temp = soln_fault(1:end,3).*pi/180;
    theta0_fault_temp(type==1)=[];
    theta0_fault{i}=theta0_fault_temp;
    
    mu_fault{i}=(theta0_fault{i}-theta0)';    
    
    from = line(loseline(i),1);
    to = line(loseline(i),2);
    xij = line(loseline(i),4);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
    fij(1) = []; %remove slack bus
    dJ = -1/xij*(fij*fij');
    Ji = J0 + dJ; %H matrix in paper
    
    if rank(Ji) < size((Ji),1)
        disp(['bad line from=' num2str(from), ', to=', num2str(to)]);
        disp('singular Jacobian')
        klcov_lines(i,1) = NaN; 
        klcov_lines(i,2) = from; 
        klcov_lines(i,3) = to; 
        continue
    end
    
    Bi = inv(Ji); %Compute each M corresponding to each line outage in A1
    
    B{i} = Bi;
    sigma{i} = Ji\(2*Psig^2*eye(nbus-1))/(Ji'); %make it symmetric to avoid numerical rounding issues
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';

%     
% %     mu{i} = J0\((theta0(from)-theta0(to))*pi/180*(1/xij)*fij);
% %     
%     klmean(i) = 1/2*(mu{i}'*inv(sigma0)*mu{i});
    klcov_lines(i,1) = 1/2 *(trace(sigma0\sigma{i}) - ...
        (nbus-1) + log(det(sigma0/sigma{i}))); 
    klcov_lines(i,2) = from; 
    klcov_lines(i,3) = to; 
end
   
%Lowest KL divergence dictates the misdetection error. We want KL divergences to be big
KL_sorted_A1=sortrows(klcov_lines) %sort KL divergnce from min to max as well as their associated lines
% disp(['min KL divergence for Area 1 = ' num2str(KL_sorted_A1(1,1))]);


%%
%Create Sample paths for line loss

linefault=1; 
% ADD=0;
% FA=Inf;
% A=100*ones(1, numel(loseline)); %Threshold

%%
%gamma=10000;
%A=nline*gamma*ones(1,numel(loseline));


faulttime=0; %use importance sampling, we want fault time =0 so that we only sample from meanshift shift distribution
num_samplepaths=1;
num_samples=50000;

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

