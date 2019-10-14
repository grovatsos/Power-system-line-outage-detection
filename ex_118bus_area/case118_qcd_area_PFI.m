clear all
close all

define_constants;
mpc = case118;
mpopt = mpoption('PF_ALG', 1, 'VERBOSE', 0, 'OUT_ALL', 0);
% [baseMVA, bus, gen, line, success, et, ybus, Yf, Yt, J] = 

results = runpf(mpc, mpopt);
baseMVA = results.baseMVA;
bus = results.bus;
line = results.branch;
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, line);

nbus = size(bus,1);
nline = size(line,1);
type = round(bus(:,2));

Ybus(type==3,:) = [];
Ybus(:,type==3) = [];

A1 = [1:45, 113:115, 117];
A2 = [46:112, 116, 118];

interarea = ismember(line(:,1),A1)&ismember(line(:,2),A2) | ...
     ismember(line(:,1),A2)&ismember(line(:,2),A1);
interarea = find(interarea==1);

% Psig = diag(0.03*ones(nbus-1,1));
% Psig = diag(0.001*bus(type~=1,4)) + ...
%         diag(0.001*bus(type~=1,6)) + ...
%         diag(0.003*ones(nbus-1,1));
Psig = diag(0.03*ones(nbus,1));

[u,I,J] = unique(line(:,1:2), 'rows');
hasDuplicates = size(u,1) < size(line(:,1:2),1);
ixDupRows = setdiff(1:size(line(:,1:2),1), I);
dupRowValues = line(ixDupRows,1:2);

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
for i = ixDupRows
    from = line(i,1);
    to = line(i,2);
    Xl = line(i,4);
    J0(from,to) = J0(from,to) - 1/Xl;
    J0(to,from) = J0(to,from) - 1/Xl;
end

%% Area 2 slack
% J11 = J0(A1,A1);
% J12 = J0(A1,A2);
% J12(:,A2==find(type==3)) = [];
% J21 = J0(A2,A1);    
% J21(A2==find(type==3),:) = [];
% J22 = J0(A2,A2);
% J22(A2==find(type==3),:) = [];
% J22(:, A2==find(type==3)) = [];

% A1_ = A1;
% A1_(A1_>(find(type==3))) = A1_(A1_>(find(type==3)))-1;
% A2_ = A2;
% A2_(A2_==find(type==3)) = [];
% A2_(A2_>(find(type==3))) = A2_(A2_>(find(type==3)))-1;

%% Area 1 slack
J11 = J0(A1,A1);
J11(A1==find(type==3),:) = [];
J11(:, A1==find(type==3)) = [];
J12 = J0(A1,A2);
J12(A1==find(type==3),:) = [];
J22 = J0(A2,A2);

A1_ = A1;
A1_(A1_==find(type==3)) = [];
A1_(A1_>(find(type==3))) = A1_(A1_>(find(type==3)))-1;
A2_ = A2;
A2_(A2_>(find(type==3))) = A2_(A2_>(find(type==3)))-1;

iflow = results.branch(interarea,PF)/baseMVA;
interflow = zeros(nbus,1);
for i = 1:length(iflow)
    interflow(line(interarea(i),1)) = interflow(line(interarea(i),1)) +...
        iflow(i);
end
% interflow = interflow(A1_);   %Area 2 slack
interflow = interflow(A1);  % Area 1 slack
interflow(type==3) = [];    % Area 1 slack

J12_tilde = sum(J12,2);
J12_tilde = diag(J12_tilde);

J0(:,type==3) = [];
J0(type==3,:) = [];

%% Compute the variance matrices corresponding to all 
%  possible line outages

B0 = inv(J11+J12_tilde);
mu0 = zeros(length(A1_),1);
sigma0 = B0*(2*Psig(1:length(A1_),1:length(A1_)).^2*eye(length(A1_)))*B0';

loseline = 1:nline;
temp = find(ismember(line(:,1),A2)|ismember(line(:,2),A2));
loseline(unique([7 9 30 59 60 61 184 temp'])) = [];
B = cell(numel(loseline),1);
sigma = cell(numel(loseline),1);
mu = cell(numel(loseline),1);
klmean = zeros(numel(loseline),1);
klcov = zeros(numel(loseline),1);

for i = 1:numel(loseline)
    
    from = line(loseline(i),1);
    to = line(loseline(i),2);
    xij = line(loseline(i),4);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
    fij(type==3) = [];    % Area 1 slack
    fij = fij(A1_);
    dJ = -1/xij*(fij*fij');
    Ji = J11 + dJ + J12_tilde;
    Bi = inv(Ji);

    B{i} = Bi;
    sigma{i} = Ji\(2*Psig(1:length(A1_),1:length(A1_)).^2*...
        eye(length(A1_)))/(Ji');
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';

    klcov(i) = 1/2 *(trace(sigma0\sigma{i}) - ...
        (length(A1_)-1) + log(det(sigma0/sigma{i})));
    
end

%% Compute relevant KL divergence values

% out = 12; % min KL div
% out = 30;
% out = 14;
% out = 7;  % max KL div
out = [12 30 14 7];
beta = [1 6 12 24 48 168];
Lbeta = length(loseline)*30*3600.*beta;
Avec = log(Lbeta);

mu11 = zeros(length(out),1);

for i = 1:length(out)
    
    outi = out(i);
    mu11(i) = computeKL(sigma{outi},sigma0,length(A1_));
    
    mu12 = zeros(length(loseline),1);
    alpha12 = zeros(length(loseline),length(Avec));
    for j = 1:length(loseline)
        mu12(j) = mu11(i) - computeKL(sigma{outi},sigma{j},length(A1_));
        alpha12(j,:) = exp(-Avec*(1 - mu12(j)/mu11(i)));
    end
    alpha12(outi,:) = zeros(1,length(Avec));
    sum(alpha12)
end

mu12 = zeros(length(loseline),length(out));
for i = 1:length(out)
    outi = out(i);
    % Compute D(fi || f0)
    mu11(i) = computeKL(sigma{outi},sigma0,length(A1_));
    for j = 1:length(loseline)
        % Compute D(fi || f0) - D(fi || fj)
        mu12(j,i) = mu11(i) - computeKL(sigma{outi},sigma{j},length(A1_));
    end
end

    