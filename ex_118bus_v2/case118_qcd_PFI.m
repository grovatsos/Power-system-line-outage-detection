clear all
close all

%% Initialize case
define_constants;
mpc = case118;
mpopt = mpoption('PF_ALG', 1, 'VERBOSE', 0, 'OUT_ALL', 0);
results = runpf(mpc, mpopt);
baseMVA = results.baseMVA;
bus = results.bus;
line = results.branch;
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, line);

nbus = size(bus,1);
nline = size(line,1);
type = round(bus(:,2));

[B0 B loseline] = makeB(case118);
nlose = length(loseline);

Psig = diag(0.03*ones(nbus,1));
Psig_ = Psig(type~=3,type~=3);

[sigma0 sigma] = makeCov(B0,B,sqrt(2)*Psig_);

KLvec = zeros(nlose,1);
for i = 1:nlose
    KLvec(i) = computeKL(sigma{i},sigma0,nbus-1);
end

%% Compute relevant KL divergence values

% KLmin = 0.1131, for i = 75
% KLmax = 2888 for i = 102
% KLmid1 = 10.2 for i = 91
% KLmid2 = 131.97 for i = 118
out = [75 91 118 102];
beta = [1 6 12 24 48 168];
Lbeta = length(loseline)*30*3600.*beta;
Avec = log(Lbeta);

mu11 = zeros(length(out),1);

for i = 1:length(out)
    
    outi = out(i);
    mu11(i) = computeKL(sigma{outi},sigma0,nbus-1);
    
    mu12 = zeros(length(loseline),1);
    alpha12 = zeros(length(loseline),length(Avec));
    for j = 1:length(loseline)
        mu12(j) = mu11(i) - computeKL(sigma{outi},sigma{j},nbus-1);
        alpha12(j,:) = exp(-Avec*(1 - mu12(j)/mu11(i)));
    end
    alpha12(outi,:) = zeros(1,length(Avec));
    sum(alpha12)
end
% 
% mu12 = zeros(length(loseline),length(out));
% for i = 1:length(out)
%     outi = out(i);
%     % Compute D(fi || f0)
%     mu11(i) = computeKL(sigma{outi},sigma0,length(A1_));
%     for j = 1:length(loseline)
%         % Compute D(fi || f0) - D(fi || fj)
%         mu12(j,i) = mu11(i) - computeKL(sigma{outi},sigma{j},length(A1_));
%     end
% end
% 
%     