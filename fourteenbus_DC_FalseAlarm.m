%%%fourteenbus_DC is the one that runs the 14bus example, generates sample paths of Wn from CuSum algorithm 
%7/16/2014: Automate some of the code so it can be applied to larger
%systems

%7/23/2014: run sample paths for simulations, see second block of this code
%8/5/2014: Note, going from 3 bus system to 14 bus system, the variance on P and Q has been changed,
           %Some lines are better for simulating outages as they affect theta significantly
           
% This code allows for PMUs at a subset of buses in addition to computing sample paths

clear all
close all
[bus, line] = buslinedata(14);

nbus = size(bus,1);
nline = size(line,1);
A1 = sort([1:nbus]); %internal buses
[soln, ybus, J] = get_pf(bus,line);
V0 = soln(1:end,2);
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
J0(:,1) = [];   J0(1,:) = []; %remove slack bus


Psig = 0.05;

B0 = inv(J0);
mu0 = zeros(length(A1)-1,1); %don't include slack
sigma0 = B0*(2*Psig^2*eye(length(A1)-1))*B0'; %M*Sigma*M' in paper, this is all internal network

loseline = 1:nline;
B = cell(numel(loseline),1);
sigma = cell(numel(loseline),1);
% mu = cell(numel(loseline),1);
klmean = zeros(numel(loseline),1);
klcov_lines = zeros(numel(loseline),3); %KL divergence and associated line loss

%Compute KL divergence for losing a line in A1
critical_line=[]; %lines where if faulted, causes islanding
for i = 1:numel(loseline)
    
    from = line(loseline(i),1);
    to = line(loseline(i),2);
    xij = line(loseline(i),4);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
    fij = fij(A1);
    fij(1) = []; %remove slack bus
    dJ = -1/xij*(fij*fij');
    Ji = J0 + dJ; %H matrix in paper
    
    if rank(Ji) < size((Ji),1)
        disp(['bad line from=' num2str(from), ', to=', num2str(to)]);
        disp('singular Jacobian A1')
        klcov_lines(i,1) = NaN; 
        klcov_lines(i,2) = from; 
        klcov_lines(i,3) = to; 
        critical_line=[critical_line i];
        continue
    end
    
    Bi = inv(Ji); %Compute each M corresponding to each line outage in A1
    
    B{i} = Bi;
    sigma{i} = Ji\(2*Psig^2*eye(length(A1)-1))/(Ji'); %make it symmetric to avoid numerical rounding issues
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';

%     
% %     mu{i} = J0\((theta0(from)-theta0(to))*pi/180*(1/xij)*fij);
% %     
%     klmean(i) = 1/2*(mu{i}'*inv(sigma0)*mu{i});
    klcov_lines(i,1) = 1/2 *(trace(sigma0\sigma{i}) - ...
        (length(A1)-1) + log(det(sigma0/sigma{i}))); 
    klcov_lines(i,2) = from; 
    klcov_lines(i,3) = to; 
end
   
%Lowest KL divergence dictates the misdetection error. We want KL divergences to be big
KL_sorted_A1=sortrows(klcov_lines); %sort KL divergnce from min to max as well as their associated lines
disp(['min KL divergence for Area 1 = ' num2str(KL_sorted_A1(1,1))]);

% %compute the D(f5||fi), want this to be large as possible to have lower false isolation rates
% %line 1 has largest D(f5||f1), thus Wn(1) for fault one line 5 will remain close to 0
% %line 2 has smallest D(f5||f2), thus Wn(2) will grow fastest for a line 5 fault
% D_f5_fi=zeros(numel(loseline),2);
% D_fi_f5=zeros(numel(loseline),2);
% for i=1:numel(loseline)
%     if ismember(i, critical_line)
%         D_f5_fi(i,1)=NaN;
%         D_f5_fi(i,2)=NaN;
%         D_fi_f5(i,1)=NaN;
%         D_fi_f5(i,2)=NaN;
%         continue;
%     end
%     if i==5
%         D_f5_fi(i,1)=0;
%         D_f5_fi(i,2)=i;
%         D_fi_f5(i,1)=0;
%         D_fi_f5(i,2)=i;
%         continue
%     end
%     D_f5_fi(i,1)=KL_compute(sigma{5}, sigma{i}); 
%     D_f5_fi(i,2)=i;
%     D_fi_f5(i,1)=KL_compute(sigma{i}, sigma{5}); 
%     D_fi_f5(i,2)=i;
% end
% D_f5_fi_sorted=sortrows(D_f5_fi);
% D_fi_f5_sorted=sortrows(D_fi_f5);
% 
% %E510 is negative, but E520 is positive
% E510=KL_compute(sigma{5}, sigma0)-KL_compute(sigma{5}, sigma{1}); %given a line 5 fault, compute the E(log(f1/f0)), this tells us if Wn will increase or remain at 0, depending on sign 
% E520=KL_compute(sigma{5}, sigma0)-KL_compute(sigma{5}, sigma{2});

%% Create sample paths for a particular line loss
seed=3;
linefault=datasample(loseline,1); %choose a random line to fault
linefault=5;
faulttime=5;%fault time
from_fault = line(loseline(linefault),1); %from bus of faulted line
to_fault = line(loseline(linefault),2);
num_samplepaths=1;
ADD=0;
% A=400; %Threshold
Beta=30;
% A=log(numel(loseline)*Beta) %Threshold
A=100;


PMU_bus=setdiff(A1, 1); %bus that has the PMU

for i = 1:numel(loseline)
    if ismember(i, critical_line)
        continue;
    end
    Bi=B{i};
    sigma{i} = Bi*(2*Psig^2*eye(length(A1)-1))*Bi'; %make it symmetric to avoid numerical rounding issues
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';
end

mu0=zeros(1, numel(theta0)); %mean of theta2, theta3....
mu = cell(numel(loseline),1);
theta0_fault=cell(numel(loseline,1));
for i=1:numel(loseline)
    if ismember(i, critical_line)
        continue;
    end
    line2=line;
    line2(i,:)=[]; %remove faulted line
    [soln_fault, ybus_fault, J_fault] = get_pf(bus,line2);
    theta0_fault{i} = soln_fault(1:end,3).*pi/180;
    from = line(loseline(i),1); %Jacobian will be singular for islanding cases, just ignore these line faults
    to = line(loseline(i),2);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
%     fij = fij(A1);
    fij(1) = []; %remove slack bus
    if to==1 %slack bus
        Pmn=(1/line(i,4))*(theta0_fault{i}(from));
    elseif from==1
        Pmn=(1/line(i,4))*(0-theta0_fault{i}(to));
    else
        Pmn=(1/line(i,4))*(theta0_fault{i}(from)-theta0_fault{i}(to));
    end
    mu{i}=-J0\(Pmn*fij); %see equations 24-25 in paper, the mean shift mu
end

for pathId=0:(num_samplepaths-1) %pathID start with 0 due to cumulative averaging formula later on
    
    tao=1; 
    f_stop=0;
    y=theta0;
    Wn=zeros(1, numel(loseline)); %CuSum, one for each line outage
    Wmun=zeros(1, numel(loseline)); %CuSum algorithm taking into account mean shift
    dyvec = [];
    Wnvec = [];
    Wmunvec = [];
    
    %initialize first instant of y. If the simulation starts off with
    %fault, then ynom is faulted y value. 
    if tao >= faulttime 
        if linefault>0 && linefault <=numel(loseline)
            y=theta0_fault{linefault};
        else
            error('invalid fault');
        end
    end
    
    
    while f_stop==0
        
        %This statement will be removed and the X vector must come from
        %Christine's programs. May be a function call?
%         dy=mvnrnd(mu, sigma0, 1);
        

%seed the random number generator
        P=-bus(2:end,6); %Negative signs are to match earlier version of code. Note that Pload and Qload may become negative meaning bus is generating power. 
        Q=-bus(2:end,7);
        P_sample=P+randn(numel(P),1)*0.05;
        Q_sample=Q+randn(numel(Q),1)*0.02;
         
        if tao < faulttime %no fault
            bus_sample=bus;
            bus_sample(2:end,6)=-P_sample;
            bus_sample(2:end,7)=-Q_sample;
            [soln_sample, ybus_sample, J_sample] = get_pf(bus_sample,line);
%             V0 = soln(1:end,2);
            theta_sample = soln_sample(1:end,3).*pi/180; %in radians
            ynew=theta_sample;
            dy = ynew'-y';
            dy = dy(2:end);
            y = ynew;
        else
            bus_sample=bus;
            bus_sample(2:end,6)=-P_sample;
            bus_sample(2:end,7)=-Q_sample;
            line2=line;
            line2(linefault,:)=[];
            [soln_sample, ybus_sample, J_sample] = get_pf(bus_sample,line2);
%             V0 = soln(1:end,2);
            theta_sample = soln_sample(1:end,3).*pi/180; %in radians
            ynew=theta_sample;
            dy = ynew'-y';
            dy = dy(2:end);
            y = ynew;
           
        end

        dyvec = [dyvec; dy];
%         X=mvnrnd(mu, sigma0, 1);

        
        %Compute the CUSUM statistics for each possible post-change
        %scenario
        
        % Algorithm 1
        for i=1:numel(loseline)
        
            if ismember(i, critical_line)
                continue;
            end
            
        Wn(i) = subplus(Wn(i) + log(mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma{i})/mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma0)));
                
        % Algorithm 2, takes into account instantaneous change
        Wmun(i) = max(subplus(Wmun(i) + ...
            log(mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma{i})/mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma0))), ...
            log(mvnpdf(dy(PMU_bus-1),mu{i}(PMU_bus-1)',sigma0)/mvnpdf(dy(PMU_bus-1), mu0(PMU_bus-1), sigma0)));
        
        end
        
        %reflect them at zero
        Wnvec = [Wnvec; Wn];
        Wmunvec = [Wmunvec; Wmun];
        
        %Check if any one of them is above A
        if sum(subplus(Wn-ones(1, numel(loseline))*A))>0
            %if yes, update the average detection delay ADD
            if pathId==0 
                ADD = tao-1;
            else
                ADD = (pathId/(pathId+1))*ADD + (tao-1)/(pathId+1); %running average to compute ADD over many sample paths
            end
            
            f_stop=1;
            
            break;
        end
        
        %Update the stopping variable till we actually stop
        tao=tao+1;
        if mod(tao, 1000)==0
            tao
            max(max(Wnvec))
        end
%         seed=seed+1;
%         seed=randi(1000,1);
    end
    if mod(pathId, 5)==0
        pathId
        ADD
    end
end

%ADD greater than Beta? False alarm if faulttime = inf
ADD
Beta
faulttime


%Compute KL divergences and also some False Isolation parameters
% KL1=KL_compute(sigma{1}, sigma0); 
% KL2=KL_compute(sigma{2}, sigma0); 
% KL3=KL_compute(sigma{3}, sigma0);
% KL_fault=KL_compute(sigma{linefault}, sigma0);
% E120=1/2*(trace(sigma0\sigma{1})-trace(sigma{2}\sigma{1})+ log(det(sigma0/sigma{2}))); 
% E130=1/2*(trace(sigma0\sigma{1})-trace(sigma{3}\sigma{1})+ log(det(sigma0/sigma{3})));
% E210=KL_compute(sigma2, sigma0)-KL_compute(sigma2, sigma1); 
% E310=1/2*(trace(sigma0\sigma3)-trace(sigma1\sigma3)+ log(det(sigma0/sigma1))); 
% E320=1/2*(trace(sigma0\sigma3)-trace(sigma2\sigma3)+ log(det(sigma0/sigma2))); 


figure;
title('Xichen Code 3bus limited PMUs')
hold on, grid on;
p5=plot(1:tao, Wnvec(1:tao,linefault)','b', 'LineWidth', 3);
p2=plot(1:tao, Wnvec(1:tao,1)','g', 'LineWidth', 3);
p3=plot(1:tao, Wnvec(1:tao,2)','r', 'LineWidth', 3);
p5_mu=plot(1:tao, Wmunvec(1:tao,linefault)','--m', 'LineWidth', 2);
p2_mu=plot(1:tao, Wmunvec(1:tao,1)','--k', 'LineWidth', 2);
p3_mu=plot(1:tao, Wmunvec(1:tao,2)', '--c', 'LineWidth', 2);
legend([p5, p2, p3, p5_mu, p2_mu, p3_mu], '5','1','2', '5 meanshift','1 meanshift','2 meanshift');



