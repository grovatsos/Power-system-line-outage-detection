%%
%Partitions the power system into two subnetwork
%Computes a sorted KL divergence of each line loss in AREA 1 (slack bus
%assumed to be in area 1)
%Computes a sorted KL divergence of each line loss in AREA 2

%Update since v2:Since the slack bus could be arbitrary chosen from the buses in Area 2,
%this code finds the maximum of min KL divergence based on this choice of
%A2 slack bus

%Intial parition is computed using grPartition.m code

clear all
close all
clc

[bus line] = buslinedata(14);
% line(29,2) = 22; %This is only for the 30 bus system to match diagram
nbus = size(bus,1);
nline = size(line,1);

% %generate a random covariance matrix with only diagonal entries
% Sigma=diag(rand(nbus-1,1));


A = adjacency_matrix(line, nbus);
[p,q,r,s] = dmperm(A); %check if graph is connected
if length(r)~=2
    error('original graph not connected')
end

% partitions=2
% [ndx,Pi,cost]= grPartition(A,partitions,1);
% 
% if ndx(1)==1
%     A1 = sort(find(ndx==1)); %internal buses
% elseif ndx(1)==2
%     A1 = sort(find(ndx==2));
% else
%     error('check partition')
% end
% A2 = sort(setdiff(1:nbus, A1)); %external buses

A1 = sort([1:5 7,8]); %internal buses
A2 = sort(setdiff(1:nbus, A1)); %external buses

if (length(A1)==1 || length(A2)==1)
    error('Too few buses in either A1 or A2')
end


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
if length(r2)~=2
    error('graph partition for A2 not connected')
end



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


Psig = 0.03;

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
disp(['min KL divergence for Area 1 = ' num2str(KL_sorted_A1(1,1))]);





%Area 2 KL divergence

max_min_KL=0;
best_A2_slack=0;
for k=1:length(A2)
    A2_slack=k;
    %Note: KL divergences change depending on which bus we chose as slack (ie. we don't measure the angle of that bus)
%     A2_slack=1; %1 doesn't mean bus 1 of overall system, but rather 1st bus in A2
    J22 = J0_orig(A2,A2);
    J22(A2_slack,:) = []; %remove slack bus of A2
    J22(:,A2_slack) = [];
    J21 = J0_orig(A2,A1);    
    J21(A2_slack,:) = [];

    J21_tilde = sum(J21,2); %row sum
    J21_tilde = diag(J21_tilde);


    Psig = 0.03;

    if rank(J22+J21_tilde) < size((J22+J21_tilde),1)
        error('bad initial partion A2');       
    end
    B0 = inv(J22+J21_tilde); 
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
    
        if  rank(Ji) < size((Ji),1)
            if (i==1&&k==1)
                disp(['bad line from=' num2str(from), ', to=', num2str(to)]);
                disp('singular Jacobian A2');
            end
            klcov_lines(i,1) = NaN; 
            klcov_lines(i,2) = from; 
            klcov_lines(i,3) = to; 
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
%     disp(['min KL divergence for Area 2 = ' num2str(KL_sorted_A2(1,1))]);
    if max_min_KL<KL_sorted_A2(1,1)
        max_min_KL=KL_sorted_A2(1,1);
        best_A2_slack=A2_slack;
    end
end
disp(['Area 2 slack bus= ' num2str(best_A2_slack(1,1))]);
disp(['min KL divergence for Area 2 = ' num2str(max_min_KL)]);