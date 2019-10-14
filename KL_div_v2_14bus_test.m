%%
%Partitions the power system into two subnetwork
%Computes a sorted KL divergence of each line loss in AREA 1 (slack bus
%assumed to be in area 1)
%Computes a sorted KL divergence of each line loss in AREA 2 

%This code investigates the case where some measurements of thetas' in Area
%1 are neglected. That is, can we detect line outages with less angle
%measurements than the full set. the variable remove takes out rows from B
%matrix.

clear all
close all


[bus line] = buslinedata(14);
% line(29,2) = 22; %This is only for the 30 bus system to match diagram
nbus = size(bus,1);
nline = size(line,1);

A = adjacency_matrix(line, nbus);

% %generate a random covariance matrix with only diagonal entries
% Sigma=diag(rand(nbus-1,1));



A1 = [1:5, 7:8]; %internal buses
A2 = setdiff(1:nbus, A1); %external buses

interarea = ismember(line(:,1),A1)&ismember(line(:,2),A2) | ...
     ismember(line(:,1),A2)&ismember(line(:,2),A1);
interarea = find(interarea==1); %interarea contains the line #'s that are tie lines
borderbusA1=intersect([line(interarea,1)' line(interarea,2)'], A1); %Node in A1 that has a tie line connected to A2
borderbusA2=intersect([line(interarea,1)' line(interarea,2)'], A2);

[soln ybus J] = get_pf(bus,line);
V0 = soln(1:end,2);
theta0 = soln(1:end,3);

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
A2_slack=1; %Chose which bus in A2 as the A2 slack bus, affects KL div of A2!
J22 = J0(A2,A2);
J22(A2_slack,:) = []; %remove slack bus of A2
J22(:,A2_slack) = [];
J21 = J0(A2,A1);    
J21(A2_slack,:) = [];


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


Psig = 0.03;


B0 = inv(J11+J12_tilde);

%%%%%%
remove=[]; %remove row of B0, equivalent to not having those angles' measurements
%%%%%%
 
B0(remove,:)=[];

mu0 = zeros(length(A1)-1,1); %don't include slack
sigma0 = B0*(2*Psig^2*eye(length(A1)-1))*B0'; %M*Sigma*M' in paper, this is all internal network

loseline = 1:nline;
temp = find(ismember(line(:,1),A2)|ismember(line(:,2),A2)); %a line that has node in A2
loseline(unique([temp'])) = []; %can't lose lines 13 or 16, else Jacobian is singular, can't lose lines with a node in A2
% loseline(unique([13 16 30 temp'])) = [];
B = cell(numel(loseline),1);
sigma = cell(numel(loseline),1);
mu = cell(numel(loseline),1);
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
    
    if rank(Ji) < length(Ji)
        disp(['bad line from=' num2str(from), ', to=', num2str(to)]);
        disp('singular Jacobian A1')
        klcov_lines(i,1) = NaN; 
        klcov_lines(i,2) = NaN; 
        klcov_lines(i,3) = NaN; 
        continue
    end
    
    Bi = inv(Ji); %Compute each M corresponding to each line outage in A1
    
    Bi(remove,:)=[];
    
    B{i} = Bi;
    sigma{i} = Bi*(2*Psig^2*eye(length(A1)-1))*Bi'; 
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';

%     
% %     mu{i} = J0\((theta0(from)-theta0(to))*pi/180*(1/xij)*fij);
% %     
%     klmean(i) = 1/2*(mu{i}'*inv(sigma0)*mu{i});
    klcov_lines(i,1) = 1/2 *(trace(sigma0\sigma{i}) - ...
        (length(sigma0)) + log(det(sigma0/sigma{i}))); 
    klcov_lines(i,2) = from; 
    klcov_lines(i,3) = to; 
end
   
%Lowest KL divergence dictates the misdetection error. We want KL divergences to be big
KL_sorted_A1=sortrows(klcov_lines) %sort KL divergnce from min to max as well as their associated lines
disp(['min KL divergence for Area 1 = ' num2str(KL_sorted_A1(1,1))]);





%Area 2 KL divergence

J21_tilde = sum(J21,2); %row sum
J21_tilde = diag(J21_tilde);


Psig = 0.03;


B0 = inv(J22+J21_tilde); %singular matrix!!!
mu0 = zeros(length(A2),1); 
sigma0 = B0*(2*Psig^2*eye(length(A2)-1))*B0'; %M*Sigma*M' in paper, this is all external network

loseline = 1:nline;
temp = find(ismember(line(:,1),A1)|ismember(line(:,2),A1)); %a line that has node in A1
loseline(unique([temp'])) = []; %Can't lose lines with a node in A1
% loseline(unique([13 16 30 temp'])) = [];
B = cell(numel(loseline),1);
sigma = cell(numel(loseline),1);
mu = cell(numel(loseline),1);
klmean = zeros(numel(loseline),1);
klcov_lines = zeros(numel(loseline),3); %KL divergence and associated line loss

%Compute KL divergence for losing a line in A2
for i = 1:numel(loseline)
    
    from = line(loseline(i),1);
    to = line(loseline(i),2);
    xij = line(loseline(i),4);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
    fij = fij(A2);
    fij(A2_slack) = []; %remove slack bus
    dJ = -1/xij*(fij*fij');
    Ji = J22 + J21_tilde + dJ; %H matrix in paper
    
    if rank(Ji) < length(Ji)
        disp(['bad line from=' num2str(from), ', to=', num2str(to)]);
        disp('singular Jacobian A2');
        klcov_lines(i,1) = NaN; 
        klcov_lines(i,2) = NaN; 
        klcov_lines(i,3) = NaN; 
        continue
    end
    
    Bi = inv(Ji); %Compute each M corresponding to each line outage in A2

    B{i} = Bi;
    sigma{i} = Ji\(2*Psig^2*eye(length(A2)-1))/(Ji'); %make it symmetric to avoid numerical rounding issues
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';

%     
% %     mu{i} = J0\((theta0(from)-theta0(to))*pi/180*(1/xij)*fij);
% %     
%     klmean(i) = 1/2*(mu{i}'*inv(sigma0)*mu{i});
    klcov_lines(i,1) = 1/2 *(trace(sigma0\sigma{i}) - ...
        (length(A2)-1) + log(det(sigma0/sigma{i}))); 
    klcov_lines(i,2) = from; 
    klcov_lines(i,3) = to; 
end
    
KL_sorted_A2=sortrows(klcov_lines); %sort KL divergnce from min to max as well as their associated lines
disp(['min KL divergence for Area 2 = ' num2str(KL_sorted_A2(1,1))]);
