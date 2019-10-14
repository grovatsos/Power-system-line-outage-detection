%%%threebus_DC is the one that runs the 3bus example, and the 2 functions it 
%uses are fun2bus and fun2bus_fault, the first is normal operation, the 2nd is
%faulted operation, and you can comment in and out some stuff to get the three different line faults

%6/18/2014: Changed the way mu, ynom are computed. 
%6/23/2014: Allows PMUs at a subset of buses (only one PMU since this is 3 bus system and one bus must be slack)

clear all
close all

global X12 X23 X13 P2 P3 Q2 Q3

X12 = 0.0504;
X23 = 0.0372;
X13 = 0.0636;

%%%%%%%%%%%%%%%%
%Choose fault
Fault=1; %1 for line 12, 2 for line 23, 3 for line 13
%%%%%%%%%%%%%%%%

A=1000; %Threshold

%%%%%%%%%%
% seed=randi([1 10000], 1, 1);  %seed the random number generator for power injections
seed=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
remove_PMU_bus=[2]; %bus 1 is slack. So either have PMU on bus 2, or bus 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PMU_bus=setdiff([2 3], remove_PMU_bus); %bus that has the PMU


J = [1/X12+1/X23 -1/X23; -1/X23 1/X13+1/X23]; %H in paper
B = inv(J); %M in paper


% L12 outage
f12 = [-1; 0];
dJ12 = -1/X12*(f12*f12'); %Delta J, equation 14 in paper
dB12 = -(J\f12*f12'/J)/(-X12 + f12'/J*f12);

% L23 outage
f23 = [1; -1];
dJ23 = -1/(X23)*(f23*f23');
dB23 = -(J\f23*f23'/J)/(-X23 + f23'/J*f23); %matrix inversion lemma
%dB23_test=(1/(X23-f23'*inv(J)*f23))*(inv(J)*f23)*(f23'*inv(J)')

% L13 outage
f13 = [0; -1];
dJ13 = -1/(X13)*(f13*f13');
dB13 = -(J\f13*f13'/J)/(-X13 + f13'/J*f13);

mu=[0 0]; %mean of theta2 or theta3
Psig = 0.5;
vari = diag([2*Psig^2, 2*Psig^2]); %variance of power injection

B_tilde12=B+dB12;
B_tilde23=B+dB23;
B_tilde13=B+dB13;

%         switch remove_PMU_bus
%             case 2
                B(remove_PMU_bus-1,:)=[];
                B_tilde12(remove_PMU_bus-1,:)=[];
                B_tilde23(remove_PMU_bus-1,:)=[];
                B_tilde13(remove_PMU_bus-1,:)=[];
%             case 3
%                 B(remove_PMU_bus-1,:)=[];
%             otherwise
%                 error('invalid PMU location');
%         end

sigma0=B*vari*B'; %variance of theta, no outage
sigma1=B_tilde12*vari*B_tilde12';
sigma2=B_tilde23*vari*B_tilde23';
sigma3=B_tilde13*vari*B_tilde13';

ADD = 0; %average detection delay

P20 = -1;
P30 = -0.9;
Q20 = -0.4;
Q30 = -0.6;

P2 = P20;
P3 = P30;
Q2 = Q20;
Q3 = Q30;
y0 = [0.01 0.01 1.01 0.99]'; %initial guess of [theta2 theta3, V2 V3]
options = optimset('Display','off');
[ynom,fval] = fsolve(@fun3bus,y0,options);
y = ynom; %[theta, V]

%mean shift at the time of outage. 1/X12*(0-ynom(1)) is P12, see
%equation 24, 26. Theta1 is 0
% mu1 = J\(1/X12*(0-ynom(1))*f12);
% mu2 = J\(1/X23*(ynom(1)-ynom(2))*f23);
% mu3 = J\(1/X13*(0-ynom(2))*f13);

%%%%%%%%%%%%
[ynom12,fval] = fsolve(@fun3bus_fault12,y0,options);
[ynom23,fval] = fsolve(@fun3bus_fault23,y0,options);
[ynom13,fval] = fsolve(@fun3bus_fault13,y0,options);
%my version, mean shift at the time of outage. 1/X12*(0-ynom(1)) is P12, see
%equation 24, 26. Theta1 is 0. mu[k] is assumed to be 0 (no net flow from internal to external)
mu1 = -J\(1/X12*(0-ynom12(1))*f12);
mu2 = -J\(1/X23*(ynom23(1)-ynom23(2))*f23);
mu3 = -J\(1/X13*(0-ynom13(2))*f13);

%%%%%%%%%%%%%%%%%%

% Trying some stuff out here
% P2v = P20 + randn(1000,1)*0.5;
% P3v = P30 + randn(1000,1)*0.5;
% Q2v = Q20 + randn(1000,1)*0.2;
% Q3v = Q30 + randn(1000,1)*0.2;
% dP2 = P2(2:end) - P2(1:end-1);
% dQ2 = Q2(2:end) - Q2(1:end-1);
% yv = zeros(4,1000);
% for i = 1:1000
%     P2 = P2v(i);
%     P3 = P3v(i);
%     Q2 = Q2v(i);
%     Q3 = Q3v(i);
%     
%     [y,fval] = fsolve(@fun3bus_fault,y0,options);
%     yv(:,i) = y;
% end

%Run for num_samplepaths sample paths
tic
num_samplepaths=1;
for pathId=0:(num_samplepaths-1) %pathID start with 0 due to cumulative averaging formula later on
    tao=1; 
    faulttime=5;%fault time
    f_stop=0;
    Wn=[0 0 0]; %CuSum, one for each line outage
    Wmun=[0 0 0]; %CuSum algorithm taking into account mean shift
    dyvec = [];
    Wnvec = [];
    Wmunvec = [];
    
    %initialize first instant of y. If the simulation starts off with
    %fault, then ynom is faulted y value. 
    if tao >= faulttime
        switch Fault
            case 1
                y=ynom12;
            case 2
                y=ynom23;
            case 3
                y=ynom13;
            otherwise
                error('invalid fault');
        end
       
    end
    
    
    while f_stop==0
        
        %This statement will be removed and the X vector must come from
        %Christine's programs. May be a function call?
%         dy=mvnrnd(mu, sigma0, 1);
        

%seed the random number generator
        rng(seed)
        P2 = P20 + randn*0.5;
        P3 = P30 + randn*0.5;
        Q2 = Q20 + randn*0.2;
        Q3 = Q30 + randn*0.2;
        
        if tao < faulttime %no fault
            options = optimset('Display','off');
            [ynew,fval] = fsolve(@fun3bus,y0,options);
            dy = ynew'-y';
            dy = dy(1:2);
            y = ynew;
        else
            options = optimset('Display','off');
            switch Fault
                case 1
                    [ynew,fval] = fsolve(@fun3bus_fault12,y0,options);
                case 2
                	[ynew,fval] = fsolve(@fun3bus_fault23,y0,options);
                case 3
                	[ynew,fval] = fsolve(@fun3bus_fault13,y0,options);
                otherwise
                    error('invalid fault');
            end
            dy = ynew'-y';
            dy = dy(1:2);
            y = ynew;
        end

        dyvec = [dyvec; dy];
%         X=mvnrnd(mu, sigma0, 1);
        
        %Compute the CUSUM statistics for each possible post-change
        %scenario
        
        % Algorithm 1
        Wn(1) = subplus(Wn(1) + log(mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma1)/mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma0)));
        Wn(2) = subplus(Wn(2) + log(mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma2)/mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma0)));
        Wn(3) = subplus(Wn(3) + log(mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma3)/mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma0))); 
        
        % Algorithm 2, takes into account instantaneous change
        Wmun(1) = max(subplus(Wmun(1) + ...
            log(mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma1)/mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma0))), ...
            log(mvnpdf(dy(PMU_bus-1),mu1(PMU_bus-1)',sigma0)/mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma0)));
        
        Wmun(2) = max(subplus(Wmun(2) + ...
            log(mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma2)/mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma0))), ...
            log(mvnpdf(dy(PMU_bus-1), mu2(PMU_bus-1)',sigma0)/mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma0)));
        
        Wmun(3) = max(subplus(Wmun(3) + ...
            log(mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma3)/mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma0))), ...
            log(mvnpdf(dy(PMU_bus-1), mu3(PMU_bus-1)',sigma0)/mvnpdf(dy(PMU_bus-1), mu(PMU_bus-1), sigma0)));

        
        %reflect them at zero
        Wnvec = [Wnvec; Wn];
        Wmunvec = [Wmunvec; Wmun];
        
        %Check if any one of them is above A
        if sum(subplus(Wn-[A A A]))>0
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
        seed=seed+1;
    end
    
end
toc
ADD


%Compute KL divergences and also some False Isolation parameters
KL12=1/2 *(trace(sigma0\sigma1)-(length(sigma0)) + log(det(sigma0/sigma1)));
KL23=1/2 *(trace(sigma0\sigma2)-(length(sigma0)) + log(det(sigma0/sigma2)));
KL13=1/2 *(trace(sigma0\sigma3)-(length(sigma0)) + log(det(sigma0/sigma3)));
E120=1/2*(trace(sigma0\sigma1)-trace(sigma2\sigma1)+ log(det(sigma0/sigma2))); 
E130=1/2*(trace(sigma0\sigma1)-trace(sigma3\sigma1)+ log(det(sigma0/sigma3)));
E210=1/2*(trace(sigma0\sigma2)-trace(sigma1\sigma2)+ log(det(sigma0/sigma1))); 
E310=1/2*(trace(sigma0\sigma3)-trace(sigma1\sigma3)+ log(det(sigma0/sigma1))); 
E320=1/2*(trace(sigma0\sigma3)-trace(sigma2\sigma3)+ log(det(sigma0/sigma2))); 


figure;
title('Xichen Code 3bus limited PMUs')
hold on, grid on;
p12=plot(1:tao, Wnvec(1:tao,1)','b', 'LineWidth', 3);
p23=plot(1:tao, Wnvec(1:tao,2)','g', 'LineWidth', 3);
p13=plot(1:tao, Wnvec(1:tao,3)','r', 'LineWidth', 3);
p12_mu=plot(1:tao, Wmunvec(1:tao,1)','--m', 'LineWidth', 2);
p23_mu=plot(1:tao, Wmunvec(1:tao,2)','--k', 'LineWidth', 2);
p13_mu=plot(1:tao, Wmunvec(1:tao,3)', '--c', 'LineWidth', 2);
legend([p12, p23, p13, p12_mu, p23_mu, p13_mu], '12','23','13', '12 meanshift','23 meanshift','13 meanshift');

% figure;
% plot(dyvec(1:tao,:));

% load('Wnvec_3bus_ex.mat');
% figure('position', [20 80 1200 800]);
% axes('position', [0.12 0.16 0.82 0.76]);
% plot(Wnvec(:,3),'b--<','LineWidth',2,'MarkerSize', 10, 'MarkerEdgeColor','b');
% hold on, grid on;
% plot(Wnvec(:,1),'g--o','LineWidth',2,'MarkerSize', 10, 'MarkerEdgeColor','g');
% plot(Wnvec(:,2),'r--s','LineWidth',2,'MarkerSize', 10, 'MarkerEdgeColor','r');
% set(gca,'fontsize',34,'fontname','times new roman'); ...
%    xl=xlabel('$k$','fontsize',38,'fontname','times new roman');
%    set(xl,'Interpreter','latex');
%    yl=ylabel('$W_k^{(m,n)}$','fontsize',38,'fontname','times new roman');
%    set(yl,'Interpreter','latex');hold on
%    ll=legend('$W_k^{(1,3)}$','$W_k^{(1,2)}$','$W_k^{(2,3)}$','Location','NorthWest');
%    set(ll,'Interpreter','latex');hold on
% axis([0.01,17,0,20]);

% load('Wn_3bus_2alg.mat');
% figure('position', [20 80 1200 800]);
% axes('position', [0.12 0.16 0.82 0.76]);
% Wnvec = [zeros(1,3); Wnvec];
% Wmunvec = [zeros(1,3); Wmunvec];
% plot(0:7, Wnvec(:,3),'b--<','LineWidth',2,'MarkerSize', 10, 'MarkerEdgeColor','b');
% hold on, grid on;
% plot(0:7, Wnvec(:,1),'g--o','LineWidth',2,'MarkerSize', 10, 'MarkerEdgeColor','g');
% plot(0:7, Wnvec(:,2),'r--s','LineWidth',2,'MarkerSize', 10, 'MarkerEdgeColor','r');
% plot(0:7, Wmunvec(:,3),'b:<','LineWidth',4,'MarkerSize', 10, 'MarkerEdgeColor','b');
% plot(0:7, Wmunvec(:,1),'g:o','LineWidth',4,'MarkerSize', 10, 'MarkerEdgeColor','g');
% plot(0:7, Wmunvec(:,2),'r:s','LineWidth',4,'MarkerSize', 10, 'MarkerEdgeColor','r');
% set(gca,'fontsize',34,'fontname','times new roman'); ...
%    xl=xlabel('$k$','fontsize',38,'fontname','times new roman');
%    set(xl,'Interpreter','latex');
%    yl=ylabel('$W_k^{(m,n)}$','fontsize',38,'fontname','times new roman');
%    set(yl,'Interpreter','latex');hold on
%    ll=legend('$W_k^{(1,3)}$','$W_k^{(1,2)}$','$W_k^{(2,3)}$',...
%        '$W_k^{(1,3)^*}$','$W_k^{(1,2)^*}$','$W_k^{(2,3)^*}$',...
%        'Location','NorthWest');
%    set(ll,'Interpreter','latex');hold on
% axis([0.01,8,0,20]);