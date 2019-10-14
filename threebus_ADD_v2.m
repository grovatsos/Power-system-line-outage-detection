% This is the code for simulating the shewhart test, cusum, meanshift, covariance shift. It samples from a
% prechange and post change distribution. No importance sampling is used. 
% Simulates the avergage ADD
%8/30/2015


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
ADD=0;
% FA=Inf;
A=45*ones(1, numel(loseline)); %Threshold

%%
%gamma=10000;
%A=log(nline*gamma)*ones(1,numel(loseline));
%%

faulttime=1; %use importance sampling
num_samplepaths=1000;
% num_samples=1000;


for pathId=1:num_samplepaths 
%     seed=pathId;
%     rng(seed);
%     Y_k=[];
    tao=1; %time
    f_stop=0;
    f_single=[]; %Specifies which single line outage streams have crossed threshold 
    dyvec = [];
    Wn_Cusum=zeros(1, numel(loseline)); %CuSum, one for each line single line outage 
    Wn_Meanshift=zeros(1, numel(loseline));
    Wn_Shewhart=zeros(1, numel(loseline));
    Wn_Covshift=zeros(1, numel(loseline));
    Wnvec_Cusum = [];
    Wnvec_Meanshift = [];
    Wnvec_Shewhart = [];
    Wnvec_Covshift = [];
    y=theta0; ynew=[]; dy=[];
    %initialize first instant of y. If the simulation starts off with fault, then ynom is faulted y value. 
    if tao >= faulttime 
       y=theta0_fault{linefault};
    end    
    
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
            dy=mvnrnd(zeros(1, numel(theta0)),sigma{linefault});
           
        end
        dyvec = [dyvec; dy];

        %Compute the CUSUM statistics for each possible post-change scenario
        for i=1:numel(loseline)
%             if ismember(loseline_single(i), critical_lines_single)
%                 continue
%             end
            Wn_Cusum(i) = max(subplus(Wn_Cusum(i) + ...
            log(mvnpdf(dy, mu0, sigma{i})/mvnpdf(dy, mu0, sigma0))), ...
            log(mvnpdf(dy,mu_fault{i},sigma0)/mvnpdf(dy, mu0, sigma0))); 
            
            Wn_Meanshift(i) = subplus(log(mvnpdf(dy,mu_fault{i},sigma0)/mvnpdf(dy, mu0, sigma0)));
            Wn_Covshift(i) = subplus(log(mvnpdf(dy, mu0, sigma{i})/mvnpdf(dy, mu0, sigma0)));
            
            Wn_Shewhart(i) = max(log(mvnpdf(dy, mu0, sigma{i})/mvnpdf(dy, mu0, sigma0)),...
            log(mvnpdf(dy,mu_fault{i},sigma0)/mvnpdf(dy, mu0, sigma0)));
        end
        Wnvec_Cusum = [Wnvec_Cusum; Wn_Cusum];
        Wnvec_Meanshift = [Wnvec_Meanshift; Wn_Meanshift];
        Wnvec_Shewhart = [Wnvec_Shewhart; Wn_Shewhart];
        Wnvec_Covshift = [Wnvec_Covshift; Wn_Covshift];
       
       
      
%       Y_k=[Y_k max(Wn_Shewhart)];
%         Y_k=[Y_k max(Wn_Covshift)];
 
%         if sum(subplus(Wn_Cusum-A))>0
%       if sum(subplus(Wn_Meanshift-A))>0
      if sum(subplus(Wn_Shewhart-A))>0
%         if sum(subplus(Wn_Covshift-A))>0
         
            [max_temp max_single]=max(Wn_Shewhart-A);
            Sorted_single=sort(subplus(Wn_Shewhart-A), 'descend');
%               find(temp(1)==(Wn_single-A))
            f_single=find(subplus(Wn_Shewhart-A)); 
%             disp(['single line outage streams ', num2str(f_single), ' crossed the treshold ', num2str(A(f_single))]);
            
            if tao<faulttime || max_single~=linefault
               error('misdetection!')
            end
            
            if pathId==1 
%                 disp(['tao = ' num2str(tao)]);
                ADD = tao-faulttime+1;
            else
%                 disp(['tao = ' num2str(tao)]);
                ADD = ((pathId-1)*ADD + (tao-faulttime+1))/pathId; %running average to compute ADD over many sample paths
            end
            pathId;
            f_stop=1;
            break; %needed so tao_single doesn't increase  
        end 
        %Update the stopping variable till we actually stop
        if mod(tao, 1000)==0
            disp(['tao = ', num2str(tao)]);
        end
        
        tao=tao+1;
       
        
        if tao>10000
            f_stop=1;
            error('timed out')
        end    
        
    end
    if mod(pathId, 10)==0
            disp(['current path = ', num2str(pathId)]);
            ADD
            tao

    end
end
 
ADD

% figure;
% plot(Wnvec_Cusum, '*');
% title('cusum');

figure;
plot(Wnvec_Shewhart);
title('shewhart')

figure;
plot(Wnvec_Covshift);
title('Covariance shift only')

figure;
plot(Wnvec_Covshift(:,linefault), '*');
title('Covariance shift only')

figure;
plot(Wnvec_Meanshift(:,linefault), '*');
title('meanshift only');

figure;
plot(Wnvec_Cusum);
title('cusum');