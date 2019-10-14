%%
%Partitions the power system into two subnetwork
%Computes a sorted KL divergence of each line loss in AREA 1 (slack bus
%assumed to be in area 1)
%Computes a sorted KL divergence of each line loss in AREA 2

%Intial parition is computed using grPartition.m code
%For thirty bus system, remember to comment out line 19

%5/17/2014 Implemented the Greedy Algorithm in finding optimal placement of
%PMUs one bus at a time. The new matrix M_tilde can be modelled as C*M
%where C has one 1 in each row and each column of C sums to one. The code
%performs exhaustive search of the column position of the 1 in C for 1 PMU,
%then another, etc.

%9/30/2014 Added code to compute the change in KL divergence for each line as additional PMUs are added. The main 
%objective of this is to see if each additional PMU ALWAYS increase the KL divergence for a certain line outage
%The variable of interest is KL_delta

clear all
close all
% clc

[bus line] = buslinedata(30);
%%%%%%%%%%%
line(29,2) = 22; %This is only for the 30 bus system to match diagram
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% use grPartition code to parition system
% [ndx,Pi,cost]= grPartition(A,partitions,1); %A program that paritions the power system 
% 
% if ndx(1)==1
%     A1 = sort(find(ndx==1)); %internal buses, let bus 1 be slack, and in area 1
% elseif ndx(1)==2
%     A1 = sort(find(ndx==2));
% else
%     error('check partition')
% end
% A2 = sort(setdiff(1:nbus, A1)); %external buses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% A1 = sort([1:5 7,8]); %internal buses, 14 bus system
A1 = sort([1:9, 11:14, 16]); %internal buses, 30 bus system
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





% %Area 2 KL divergence
% %k2 is number of PMUs allocated to Area 2
% 
% %%%%%%
% %Specify number of PMUs allocated for Area 2, must be less than length of A2 since at least 1 bus must be slack
% k2=length(A2)-6; 
% %%%%%%
% 
% if k2>length(A2)-1 || k2<1
%     error('num PMUs allocated error')
% end
% 
% C=nchoosek(A2,k2); %all combinations of buses with PMU measurements
% max_min_KL=0;
% KL_sorted_A2=cell(3,length(C)*(length(A2)-k2)); %Structure of this cell is KL divergences for all removed buses combinations and then for each combination, the bus desginated as slack
% optimal_remove=[];
% count=1;
% for k=1:length(C)
%      bus_no_PMUs=setdiff(A2,C(k,:));
%      remove=find(ismember(A2,bus_no_PMUs));
%     %Note: KL divergences change depending on which bus we chose as slack (ie. we don't measure the angle of that bus)
%     
%     %iterate through all possible slack bus locations over all elements of remove
%     for j=1:length(remove)
%     
%     A2_slack=remove(j); %1 doesn't mean bus 1 of overall system, but rather 1st bus in A2
%     J22 = J0_orig(A2,A2);
%     J22(A2_slack,:) = []; %remove slack bus of A2
%     J22(:,A2_slack) = [];
%     J21 = J0_orig(A2,A1);    
%     J21(A2_slack,:) = [];
% 
%     J21_tilde = sum(J21,2); %row sum
%     J21_tilde = diag(J21_tilde);
% 
% 
%     Psig = 0.03;
% 
%     if rank(J22+J21_tilde) < size((J22+J21_tilde),1)
%         error('bad initial partion A2');       
%     end
%     B0 = inv(J22+J21_tilde); 
%    
%     if length(remove)>1
%        B0(remove(find((remove-A2_slack)>0))-1,:)=[];
%        B0(remove(find((remove-A2_slack)<0)),:)=[]; 
% %        B0(remove(2:end)-1,:)=[]; %delete rows of B0, different from slack bus, which removes row and column from J, shift of 1 is due to slack bus deleteing a row and col of J22
%    end
%     mu0 = zeros(length(A2),1); 
%     sigma0 = B0*(2*Psig^2*eye(length(A2)-1))*B0'; %M*Sigma*M' in paper, this is all external network
% 
%     loseline = 1:nline;
%     temp = find(ismember(line(:,1),A1)|ismember(line(:,2),A1)); %a line that has node in A1
%     loseline(unique([temp'])) = []; %Can't lose lines with a node in A1
%     % loseline(unique([13 16 30 temp'])) = [];
%     B = cell(numel(loseline),1);
%     sigma = cell(numel(loseline),1);
%     mu = cell(numel(loseline),1);
%     klmean = zeros(numel(loseline),1);
%     klcov_lines = zeros(numel(loseline),3); %KL divergence and associated line loss
% 
%     %Compute KL divergence for losing a line in A2
%     for i = 1:numel(loseline)
%     
%         from = line(loseline(i),1);
%         to = line(loseline(i),2);
%         xij = line(loseline(i),4);
%         fij = zeros(nbus,1);
%         fij(from) = 1;
%         fij(to) = -1;
%         fij = fij(A2);
%         fij(A2_slack) = []; %remove slack bus
%         dJ = -1/xij*(fij*fij');
%         Ji = J22 + J21_tilde + dJ; %H matrix in paper
%     
%         if  rank(Ji) < size((Ji),1)
%             if (i==1&&k==1)
%                 disp(['bad line from=' num2str(from), ', to=', num2str(to)]);
%                 disp('singular Jacobian A2');
%             end
%             klcov_lines(i,1) = NaN; 
%             klcov_lines(i,2) = from; 
%             klcov_lines(i,3) = to; 
%             continue
%         end
%     
%         Bi = inv(Ji); %Compute each M corresponding to each line outage in A2
%         
%         if length(remove)>1
%             Bi(remove(find((remove-A2_slack)>0))-1,:)=[];
%             Bi(remove(find((remove-A2_slack)<0)),:)=[];
% %             Bi(remove(2:end)-1,:)=[]; %delete rows of B0, different from slack bus, which removes row and column from J
%         end
%         
%         B{i} = Bi;
%         sigma{i} = Bi*(2*Psig^2*eye(length(A2)-1))*Bi'; 
%         %make it symmetric to avoid numerical rounding issues
%         aa = triu(sigma{i});
%         sigma{i} = aa + triu(aa,1)';
% 
% %     
% % %     mu{i} = J0\((theta0(from)-theta0(to))*pi/180*(1/xij)*fij);
% % %     
% %     klmean(i) = 1/2*(mu{i}'*inv(sigma0)*mu{i});
%         klcov_lines(i,1) = 1/2 *(trace(sigma0\sigma{i}) - ...
%             (length(sigma0)) + log(det(sigma0/sigma{i}))); 
%         klcov_lines(i,2) = from; 
%         klcov_lines(i,3) = to; 
%     end
%     
%     KL_sorted_A2{1,count}=sortrows(klcov_lines); %sort KL divergnce from min to max as well as their associated lines
%     KL_sorted_A2{2,count}=A2(A2_slack);
%     if length(remove)>1
%         KL_sorted_A2{3,count}=A2(setdiff(remove, A2_slack));
%     end
%     %     disp(['min KL divergence for Area 2 = ' num2str(KL_sorted_A2(1,1))]);
%     
% 
%     if max_min_KL<KL_sorted_A2{1,count}(1,1)
%         max_min_KL=KL_sorted_A2{1,count}(1,1);
%         optimal_slack=A2(A2_slack);
%         temp=remove;
%         temp(j)=[];
%         optimal_remove=A2(temp);
%     end
%     count=count+1;
%     end
% end
% disp(sprintf(['\nnum PMUs allocated for Area 2= ' num2str(k2)]));
% disp(sprintf(['min KL divergence optimal placement of PMUs for Area 2 = ' num2str(max_min_KL), '\nslack bus = ' num2str(optimal_slack), '\nother buses without PMUs = ', num2str(optimal_remove)]));
% disp(['A2= ' num2str(A2)]);
% 
% save('KL_div.mat')



%%%%%%%%% Greedy Algorithm
%%
%Area 2 KL divergence
%k2 is number of PMUs allocated to Area 2

%%%%%%
%Specify number of PMUs allocated for Area 2, must be less than length of A2 since at least 1 bus must be slack
k2=length(A2)-1; 
%%%%%%

tic
if k2>length(A2)-1 || k2<1
    error('num PMUs allocated error')
end

% C=zeros(k2,length(A2)-1);
max_min_KL=[];
KL_sorted_A2_greedy=cell(3,k2); %structure of this cell is list optimal KL divergences as each successive row is added (each PMU added). 
% optimal_remove=[];
optimal_PMU_indicies=[];
optimal_slack=[];
optimal_k=0;

KL_delta=cell(1,k2-1); %keep change of KL changes as each PMU is successively added

for k=1:k2 %add each PMU successively
    max_min_KL_k=0;
    optimal_remove_k=[];
    optimal_slack_k=[];
    for l=1:length(A2)         %exhaustively search next PMU location by computing KL divergence
%         if isempty(find(C(:,l),1))
%             C(k,:)=zeros(1,length(A2)-1);
%             C(k,l)=1;
%         end
%         [row col]=find(C(1:k,:));
        if ismember(l,optimal_PMU_indicies)
            continue
        end
        bus_no_PMUs=setdiff(A2,A2([optimal_PMU_indicies l]));
        remove=find(ismember(A2,bus_no_PMUs));
    %Note: KL divergences change depending on which bus we chose as slack (ie. we don't measure the angle of that bus)
    
    %iterate through all possible slack bus locations over all elements of remove
        for j=1:length(remove)
    
        A2_slack=remove(j); %1 doesn't mean bus 1 of overall system, but rather 1st bus in A2
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
   
        if length(remove)>1
            B0(remove(find((remove-A2_slack)>0))-1,:)=[];
            B0(remove(find((remove-A2_slack)<0)),:)=[]; 
%           B0(remove(2:end)-1,:)=[]; %delete rows of B0, different from slack bus, which removes row and column from J, shift of 1 is due to slack bus deleteing a row and col of J22
        end
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
        
            if length(remove)>1
                Bi(remove(find((remove-A2_slack)>0))-1,:)=[];
                Bi(remove(find((remove-A2_slack)<0)),:)=[];
%               Bi(remove(2:end)-1,:)=[]; %delete rows of B0, different from slack bus, which removes row and column from J
            end
        
            B{i} = Bi;
            sigma{i} = Bi*(2*Psig^2*eye(length(A2)-1))*Bi'; 
            %make it symmetric to avoid numerical rounding issues
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
        
        
        sorted_KL=sortrows(klcov_lines); %sort KL divergnce from min to max as well as their associated lines
%         disp(sprintf(['k= ' num2str(k)]));
%         disp(sprintf(['l= ' num2str(l)]));
%         disp(sprintf(['j= ' num2str(j)]));

        if sorted_KL(1,1)<-1e-5
            error('KL div. negative')
        end
        if max_min_KL_k<sorted_KL(1,1)
            max_min_KL_k=sorted_KL(1,1);
            temp=remove;
            temp(j)=[];
            optimal_remove_k=A2(temp);
%             optimal_PMU_bus_k=A2(l);
            optimal_slack_k=A2(A2_slack);
            KL_sorted_A2_k=sorted_KL;
            optimal_PMU_bus_index=l;
            
            KL_current=klcov_lines; %unsorted KL divergences
            
        end
        end
        
    end
    
    if k>1 %for computing the increase in individual KL divergences from adding a PMU
       KL_delta{1,k-1}=KL_current(:,1)-KL_prev(:,1);
    end
    KL_prev=KL_current; %save previous KL divergence (unsorted) of lines for KL_delta computation
    

       
%         C(k,:)=zeros(1,length(A2)-1);
%         C(k,find(A2==optimal_PMU_bus_k))=1;
        optimal_PMU_indicies=[optimal_PMU_indicies  optimal_PMU_bus_index];
        KL_sorted_A2_greedy{1,k}=KL_sorted_A2_k; %sort KL divergnce from min to max as well as their associated lines
        KL_sorted_A2_greedy{2,k}=optimal_slack_k;
        if length(remove)>1
           KL_sorted_A2_greedy{3,k}=optimal_remove_k; %Buses with no PMU
        end
        
%         if max_min_KL(k)<max_min_KL_k
            max_min_KL=[max_min_KL max_min_KL_k];
            optimal_slack=[optimal_slack optimal_slack_k];
%             optimal_remove=[optimal_remove optimal_remove_k];
            
%         end
%     if mod(l, 10)==0
%        disp(sprintf(['count=' num2str(count)]));
       disp(sprintf(['k=' num2str(k)]));
%     end 
end
toc
[temp optimal_k]=max(max_min_KL);

   
% disp(sprintf(['\nmax allocated # of PMUs for Area 2= ' num2str(k2), '. The optimal number of PMUs is ', num2str(optimal_k), '\nwith PMUs at buses  ', num2str(sort(A2(optimal_PMU_indicies(1:optimal_k))))]));
% disp(sprintf(['\nmin KL divergence optimal placement of PMUs for Area 2 = ' num2str(max_min_KL(optimal_k)), '\nslack bus = ' num2str(optimal_slack(optimal_k)), '\nother buses without PMUs = ', num2str(KL_sorted_A2_greedy{3,optimal_k})]));
% disp(sprintf(['\nideal num of PMUs for Area2 is ' num2str(optimal_k)]));
% disp(['\nC= ' num2str(C)]);
% disp(['A2= ' num2str(A2)]);


save('KL_div_greedy.mat') %To see the max_min_KL for each k2 allocation of PMUs, see KL_sorted_A2_greedy{1,k}(1,1) 



%%

% figure('position', [20 80 1300 800]);
% axes('position', [0.14 0.14 0.8 0.8]);
% hold on; box on; grid on;
% 
% x=[1:19];
% y1=[ 0.000142 0.00277 0.030007 0.1677 0.442 0.9562 1.1168 1.235 1.4303... 
%     1.5786 1.6305 1.7609 1.7609 1.7609 1.7609 1.7609 1.7609 1.7609 1.7609];
% y2=[0.000142 0.00234 0.00675 0.02346 0.1748 0.4929 1.1159 1.2211 1.2274 1.3723 1.4320 1.5083 1.516 1.522 1.522 1.522 1.522 1.522 1.522];
% p1=plot(x, y1, 'r', 'LineWidth', 2, 'Marker', '*', 'MarkerSize', 15); 
% p2=plot(x, y2, 'b', 'LineWidth', 2, 'Marker', 'X', 'MarkerSize', 15); 
% 
% 
%    
% set(gca,'fontsize',36,'fontname','times new roman'); ...
% %    xl=xlabel('$ V_{110} $ [p.u.] ','fontsize',40,'fontname','times new roman');
%    xl=xlabel('No. of PMUs','fontsize',40,'fontname','times new roman');
%    set(xl,'Interpreter','latex');
%    yl=ylabel('Min. KL Divergence','fontsize',40,'fontname','times new roman');
%    set(yl,'Interpreter','latex');
%    l1=legend([p1, p2], 'Exhaustive Search', 'Greedy Algorithm' , 'location', 'SouthEast');
%    set(l1,'FontSize',37, 'Interpreter','latex')

% axis([0.945 1.055 0.945 1.055]);
