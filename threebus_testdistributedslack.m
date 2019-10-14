%%Code used for testing out solutions of Matpower and get_pf

% used for testing distributed slack bus
% clear all
% clc
% close all


[bus, line] = buslinedata(14);
% P_sample=[1.6 0.9];
% Q_sample=[0.4 0.6];
% lineout=2;
alpha=[1 1 1 0 0 1 0 1 0 0 0 0 0 0]'/5; %Distributed slack bus
alpha=[1 0 0 0 0 0 0 0 0 0 0 0 0 0]';
alpha=[0 0 0 0 0 0 0 1 0 0 0 0 0 0]';
% alpha=[1 0 0]';


% bus(2:end,6)= P_sample;
% bus(2:end,7)= Q_sample;
% line(lineout,:)=[];

[soln, ybus, J_getpf] = get_pf_v4(bus,line);
V0_getpf = soln(1:end,2)
theta0_getpf = soln(1:end,3).*pi/180 %in radians

% [soln ybus dVdQ H K N L]=get_pf_v4(bus,line); %H is the decoupled power flow jacobian




system=case14;
mpc = system;
mpopt = mpoption('PF_ALG', 1, 'VERBOSE', 0, 'OUT_ALL', 0);
nbus = size(mpc.bus,1);
nline = size(mpc.branch,1);

% bus=mpc.bus;
% bus(2:end,3)=P_sample*100;
% bus(2:end,4)=Q_sample*100;
% mpc.bus=bus;
% 
% line=mpc.branch;
% line(lineout,:)=[];
% mpc.branch=line;


results = runpf(mpc, mpopt);
V0_Matpower = results.bus(:,8)
theta0_Matpower = results.bus(:,9).*pi/180 %in radians

%The get_pf_dsb code takes input of power injections in perunit 
temp=mpc.bus;
temp(:,3:4)=temp(:,3:4)/mpc.baseMVA;
mpc.bus=temp;
temp=mpc.gen;
temp(:,2)=temp(:,2)/mpc.baseMVA;
mpc.gen=temp;

mpc.gen(:,5) = -9999; %Qmin and Qmax
mpc.gen(:,4) = 9999;

mpc.bus(:,12) = 0.9;

temp=mpc.bus;   %scale the shunts into p.u.
temp(:,5:6)=temp(:,5:6)/mpc.baseMVA
mpc.bus=temp;

[soln, ybus, J_getpf, Jfull] = get_pf_dsb(mpc,alpha);
V0_getpf_dsb = soln(1:end,2)
theta0_getpf_dsb = soln(1:end,3).*pi/180 %in radians

[V0_getpf V0_Matpower V0_getpf_dsb]
[theta0_getpf theta0_Matpower theta0_getpf_dsb]

bus=mpc.bus;
type = round(bus(:,2));
(results.bus(:,9) - results.bus(1,9))*pi/180;

J_Matpower=makeJac(mpc);



%Compare these results with the output of file threebus_DC_v3.m where the solutions is solved exactly


