%%
%Partitions the power system into two subnetwork
%Computes a sorted KL divergence of each line loss in AREA 1 (slack bus
%assumed to be in area 1)

clear all
close all

[bus line] = buslinedata(30);
line(29,2) = 22; %This is only for the 30 bus system to match diagram
nbus = size(bus,1);
nline = size(line,1);

A = incidence_matrix(line, nbus);

%generate a random covariance matrix with only diagonal entries
Sigma=diag(rand(nbus-1,1));



A1 = [1:7,9:24]; %internal buses
A2 = setdiff(1:nbus, A1); %external buses
%A2 = [8,25:30]; %external buses
interarea = ismember(line(:,1),A1)&ismember(line(:,2),A2) | ...
     ismember(line(:,1),A2)&ismember(line(:,2),A1);
interarea = find(interarea==1); %interarea contains the line #'s that are tie lines

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

J11 = J0(A1,A1);    J11(:,1) = [];  J11(1,:) = []; %bus 1 is assumed slack bus
J12 = J0(A1,A2);    J12(1,:) = [];
J21 = J0(A2,A1);    J21(:,1) = [];
J22 = J0(A2,A2);

J0(:,1) = [];   J0(1,:) = []; %remove slack bus

S = get_line_flow(line, soln, ybus);
iflow = real(S(interarea));
interflow = zeros(nbus,1); %line flow assigned to "from" bus
for i = 1:length(iflow)
    interflow(line(interarea(i),1)) = interflow(line(interarea(i),1)) +...
        iflow(i);
end
interflow = interflow(A1); %internal area tie lines flows, assigned to from bus

J12_tilde = sum(J12,2); %row sum
J12_tilde = diag(J12_tilde);


Psig = 0.03;


B0 = inv(J11+J12_tilde);
mu0 = zeros(length(A1)-1,1); %don't include slack
sigma0 = B0*(2*Psig^2*eye(length(A1)-1))*B0'; %M*Sigma*M' in paper, this is all internal network

loseline = 1:nline;
temp = find(ismember(line(:,1),A2)|ismember(line(:,2),A2)); %a line that has node in A2
loseline(unique([13 16 temp'])) = []; %can't lose lines 13 or 16, else Jacobian is singular, can't lose lines with a node in A2
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
%     if rank(Ji) < length(Ji)
%         i
%     end
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
    
KL_sorted=sortrows(klcov_lines); %sort KL divergnce from min to max as well as their associated lines
disp(['min KL divergence for Area 1 = ' num2str(KL_sorted(1,1))]);
