%%%fourteenbus_DC is the one that runs the 14bus example, generates sample paths of Wn from CuSum algorithm 
%7/16/2014: Automate some of the code so it can be applied to larger
%systems

%7/23/2014: run sample paths for simulations, see second block of this code
%8/5/2014: Note, going from 3 bus system to 14 bus system, the variance on P and Q has been changed,
           %Some lines are better for simulating outages as they affect theta significantly
           
%10/28/2014: This code checks ADD by setting individual thresholds for each line outage scaled according to their KL divergence

%11/18/2014: Adapted from 3bus to 14bus system. 
           
% This code allows for PMUs at a subset of buses in addition to computing sample paths

clear all
% close all


[bus, line] = buslinedata(14);
%%%%%%%%%%%
% line(29,2) = 22; %This is only for the 30 bus system to match diagram
%%%%%%%%%%%

nbus = size(bus,1);
nline = size(line,1);

% %generate a random covariance matrix with only diagonal entries
% Sigma=diag(rand(nbus-1,1));


A = adjacency_matrix(line, nbus);
[p,q,r,s] = dmperm(A); %check if graph is connected
if length(r)~=2
    error('original graph not connected')
end

partitions=2;

A1 = sort([1:nbus]); %internal buses
A2 = sort(setdiff(1:nbus, A1)); %external buses

% if (length(A1)==1 || length(A2)==1)
%     error('Too few buses in either A1 or A2')
% end


interarea = ismember(line(:,1),A1)&ismember(line(:,2),A2) | ...
     ismember(line(:,1),A2)&ismember(line(:,2),A1);
interarea = find(interarea==1); %interarea contains the line #'s that are tie lines

area1 = ismember(line(:,1),A1)&ismember(line(:,2),A1); %all lines in Area 1
area2 = ismember(line(:,1),A2)&ismember(line(:,2),A2); 
area1 = find(area1==1);
area2 = find(area2==1);

borderbusA1=intersect([line(interarea,1)' line(interarea,2)'], A1); %Node in A1 that has a tie line connected to A2
borderbusA2=intersect([line(interarea,1)' line(interarea,2)'], A2);


%check if partitions A1 and A2 are connected graphs
lineA1=line(area1,:);  %lines in Area 1
lineA2=line(area2,:);  %lines in Area 2

%renumber the buses of lineA1 and lineA2 to match that of A1_adj matrix
%this is needed to make adjacency function work
diffA1=A1-(1:length(A1));
indexA1=find(diffA1);
lineA1_from=lineA1(:,1)'; 
lineA1_to=lineA1(:,2)';
for i=1:length(indexA1)
    lineA1_from(lineA1_from==A1(indexA1(i)))=indexA1(i);
    lineA1_to(lineA1_to==A1(indexA1(i)))=indexA1(i);
end
lineA1_renum=[lineA1_from; lineA1_to]';

diffA2=A2-(1:length(A2));
indexA2=find(diffA2);
lineA2_from=lineA2(:,1)'; 
lineA2_to=lineA2(:,2)';
for i=1:length(indexA2)
    lineA2_from(lineA2_from==A2(indexA2(i)))=indexA2(i);
    lineA2_to(lineA2_to==A2(indexA2(i)))=indexA2(i);
end
lineA2_renum=[lineA2_from; lineA2_to]';

%compute adjacency matrices for two areas
A1_adj = adjacency_matrix(lineA1_renum, length(A1));
A2_adj = adjacency_matrix(lineA2_renum, length(A2));

[p1,q1,r1,s1] = dmperm(A1_adj); %check if graph of A1 connected
[p2,q2,r2,s2] = dmperm(A2_adj); %check if graph of A2 connected
if length(r1)~=2
    error('graph partition for A1 not connected')
end
% if length(r2)~=2
%     error('graph partition for A2 not connected')
% end

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

%internal 
J11 = J0(A1,A1);    
J11(:,1) = []; %remove slack bus equations
J11(1,:) = []; %bus 1 is assumed slack bus
J12 = J0(A1,A2);    
J12(1,:) = [];

%external 
J0_orig=J0;
J0(:,1) = [];   J0(1,:) = []; %remove slack bus

S = get_line_flow(line, soln, ybus);
iflow = real(S(interarea));
interflow = zeros(nbus,1); %line flow assigned to "from" bus
for i = 1:length(iflow)
    interflow(line(interarea(i),1)) = interflow(line(interarea(i),1)) +...
        iflow(i);
end
interflow = interflow(A1); %internal area tie lines flows, assigned to from bus

%Area 1 KL divergence
J12_tilde = sum(J12,2); %row sum
J12_tilde = diag(J12_tilde);


Psig = 0.05; %use 0.05 for 14 bus, 0.1 for 3 bus, this does affect KL divergence so in theory, it shouldn't affect ADD
%but, if it hits some limit in power flow and goes to a different operating point the H matrices may change completely 
Qsig=0.02;

if rank(J11+J12_tilde) < size((J11+J12_tilde),1)
       error('bad initial partion A1');       
end
B0 = inv(J11+J12_tilde);
mu0 = zeros(length(A1)-1,1); %don't include slack
sigma0 = B0*(2*Psig^2*eye(length(A1)-1))*B0'; %M*Sigma*M' in paper, this is all internal network

loseline = 1:nline;
temp = find(ismember(line(:,1),A2)|ismember(line(:,2),A2)); %a line that has node in A2
loseline(unique([temp'])) = []; %can't lose lines 13 or 16, else Jacobian is singular, can't lose lines with a node in A2
% loseline(unique([13 16 30 temp'])) = [];
B = cell(numel(loseline),1);
sigma = cell(numel(loseline),1);
% mu = cell(numel(loseline),1);
klmean = zeros(numel(loseline),1);
klcov_lines = zeros(numel(loseline),3); %KL divergence and associated line loss

%Compute KL divergence for losing a line in A1
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
    Ji = J11 + J12_tilde + dJ; %H matrix in paper
    
    if rank(Ji) < size((Ji),1)
        disp(['bad line from=' num2str(from), ', to=', num2str(to)]);
        disp('singular Jacobian A1')
        klcov_lines(i,1) = NaN; 
        klcov_lines(i,2) = from; 
        klcov_lines(i,3) = to; 
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
%disp(['min KL divergence for Area 1 = ' num2str(KL_sorted_A1(1,1))]);

critical_lines=find(isnan(klcov_lines)); %lines that if lost, creating islanding situations

%% Create sample paths for a particular line loss
% seed=4;

%close all
linefault=datasample(loseline,1); %choose a random line to fault
linefault=6;
if ismember(linefault, critical_lines)
    error('choose different line to fault')
end
from_fault = line(loseline(linefault),1); %from bus of faulted line
to_fault = line(loseline(linefault),2);
num_samplepaths=1;
y=theta0;
ADD=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%B0 and sigma0 are modified according to buses without PMUs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
remove_PMU_bus=[]; %buses that do not have PMUs, bus 1 is slack, so remove PMU can only remove from 2 or 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute the KL divergences with certain bus' PMUs removed

PMU_bus=setdiff([A1], [1 remove_PMU_bus]); %bus that has the PMU
B0(remove_PMU_bus-1,:)=[];
sigma0 = B0*(2*Psig^2*eye(length(A1)-1))*B0';
klcov_lines_limitedPMUs = zeros(numel(loseline),3); %KL divergence and associated line loss with less than full set of PMUs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the new KL divergences with limited number of PMUs
for i = 1:numel(loseline)
    from = line(loseline(i),1);
    to = line(loseline(i),2);
    klcov_lines_limitedPMUs(i,2) = from; 
    klcov_lines_limitedPMUs(i,3) = to; 
    if ismember(i, critical_lines)
         klcov_lines_limitedPMUs(i,1) =NaN;
        continue
    end
    Bi=B{i};
    Bi(remove_PMU_bus-1,:)=[];
    B{i} = Bi;
    sigma{i} = Bi*(2*Psig^2*eye(length(A1)-1))*Bi'; %make it symmetric to avoid numerical rounding issues
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';
    klcov_lines_limitedPMUs(i,1) = KL_compute(sigma{i}, sigma0);   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scale thresholds according to KL divergences
%run sample paths
A=(100*klcov_lines_limitedPMUs(:,1)/max(klcov_lines_limitedPMUs(:,1)))';
A=100*ones(1, numel(loseline)); %Threshold


mu0=zeros(1, numel(theta0)); %mean of theta2, theta3....
mu = cell(numel(loseline),1);
theta0_fault=cell(numel(loseline,1));

for i=1:numel(loseline)
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
    faulttime=1;%fault time
    f_stop=0;
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
%         rng(seed)
        P=-bus(2:end,6); %Negative signs are to match earlier version of code. Note that Pload and Qload may become negative meaning bus is generating power. 
        Q=-bus(2:end,7);
        P_sample=P+randn(numel(P),1)*Psig;
        Q_sample=Q+randn(numel(Q),1)*Qsig;
         
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
            line2([2 linefault],:)=[];
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
            if ismember(i, critical_lines)
                continue
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
        if sum(subplus(Wn-A))>0
            %if yes, update the average detection delay ADD
            if pathId==0 
                ADD = tao-1;
            else
                ADD = (pathId/(pathId+1))*ADD + (tao-1)/(pathId+1) %running average to compute ADD over many sample paths
            end
            
            f_stop=1;
            
            break;
        end
        
        %Update the stopping variable till we actually stop
        tao=tao+1;
%         seed=seed+1;
    end
    
end
ADD=ADD-faulttime
% disp(ADD)


%Compute KL divergences and also some False Isolation parameters
KL1=KL_compute(sigma{1}, sigma0); 
KL2=KL_compute(sigma{2}, sigma0); 
KL3=KL_compute(sigma{3}, sigma0);
KL_fault=KL_compute(sigma{linefault}, sigma0);
E120=1/2*(trace(sigma0\sigma{1})-trace(sigma{2}\sigma{1})+ log(det(sigma0/sigma{2}))); 
E130=1/2*(trace(sigma0\sigma{1})-trace(sigma{3}\sigma{1})+ log(det(sigma0/sigma{3})));
E210=KL_compute(sigma{2}, sigma0)-KL_compute(sigma{2}, sigma{1}); 
E230=KL_compute(sigma{2}, sigma0)-KL_compute(sigma{2}, sigma{3}); 
E310=1/2*(trace(sigma0\sigma{3})-trace(sigma{1}\sigma{3})+ log(det(sigma0/sigma{1}))); 
E320=1/2*(trace(sigma0\sigma{3})-trace(sigma{2}\sigma{3})+ log(det(sigma0/sigma{2}))); 


figure('position', [20 80 1300 800]);
axes('position', [0.14 0.14 0.8 0.8]);
hold on; box on; grid on;
% title('Xichen Code 3bus limited PMUs')
% hold on, grid on;
p1=plot(1:tao, Wnvec(1:tao,linefault)','--b', 'LineWidth', 2, 'Marker', '*', 'Markersize', 6);
p2=plot(1:tao, Wnvec(1:tao,2)','--g', 'LineWidth', 2,'Marker', 'o','Markersize', 6);
p3=plot(1:tao, Wnvec(1:tao,3)','--r', 'LineWidth', 2,'Marker', 's','Markersize', 6);
% p1_mu=plot(1:tao, Wmunvec(1:tao,linefault)','--m', 'LineWidth', 2);
% p2_mu=plot(1:tao, Wmunvec(1:tao,2)','--k', 'LineWidth', 2);
% p3_mu=plot(1:tao, Wmunvec(1:tao,3)', '--c', 'LineWidth', 2);
l1=legend([p1, p2, p3], '$W_{linefault}$','$W_{(1,3)}$','$W_{(2,3)}$');

set(gca,'fontsize',36,'fontname','times new roman'); ...
   xl=xlabel('$k$','fontsize',38,'fontname','times new roman');
   set(xl,'Interpreter','latex');
   yl=ylabel('$W[k]$','fontsize',38,'fontname','times new roman');
   set(yl,'Interpreter','latex');
   set(l1,'Interpreter','latex');
%    l1=legend([p1, p2, p3, p4], 'Minimum Vol.', 'Exact Solution', 'Linear Approx.', 'Voltage Constraints', 'location', 'SouthEast');
   set(l1,'FontSize',37, 'Interpreter','latex')
%    set(gca,'XTick',1.0487:0.0001:1.0492)
%    set(gca,'YTick',4952270:40:4952290)

   % axis([1.04901 1.0504 0.98 1.07])
