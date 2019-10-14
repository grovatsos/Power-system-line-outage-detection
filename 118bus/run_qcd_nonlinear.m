%run sample paths using Matpower code

% clear all
close all

%% Initialize case
define_constants;
mpc = case14;
mpopt = mpoption('PF_ALG', 1, 'VERBOSE', 0, 'OUT_ALL', 0);
results = runpf(mpc, mpopt);
baseMVA = results.baseMVA;
bus = results.bus;
line = results.branch;
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, line);

nbus = size(bus,1);
nline = size(line,1);
type = round(bus(:,2));

[J0 B0 B loseline] = makeB(mpc);
nlose = length(loseline);

Psig = diag(0.03*ones(nbus,1));
Psig_ = Psig(type~=3,type~=3);

[sigma0 sigma] = makeCov(B0,B,sqrt(2)*Psig_);

KLvec = zeros(nlose,1);
for i = 1:nlose
    KLvec(i) = KL_compute(sigma{i},sigma0);
%     function D = computeKL(post, pre, k)
% 
%     D = 1/2*(trace(pre\post) - k + log(det(pre/post)));
% 
% end
end
% KLmin = 0.1131, for i = 75
% KLmax = 2888 for i = 102
% KLmid1 = 10.2 for i = 91
% KLmid2 = 131.97 for i = 118
%%
beta = [1 6 12 24 48 168];
Lbeta = nlose*30*3600.*beta;
Avec = log(Lbeta);  % A threshold corresponding to FARs
ADDvec = zeros(1,length(Avec));
PFI_g0 = zeros(1,length(Avec));
mu0 = zeros(nbus-1,1);
nPath = 1;
isolated_g0 = zeros(nPath+1,1);
line_IDed = zeros(nPath+1, 10);

FA = 0;
out = 11; %102
line_post = line;
if ~FA
    line_post(loseline(out),:) = [];
end

kk = 1;
tic
while kk <= length(Avec)
    
    A = Avec(kk);
    Tthreshold = A*ones(1,nlose);
    
    for pathId=0:nPath

        tao=1;%stopping time
        f_stop=0;

        Wn=zeros(1,nlose);
        Wnvec = zeros(1e3,nlose);

        while f_stop == 0

            A1 = getMeas(mpc,Psig,line_post);
            A2 = getMeas(mpc,Psig,line_post);
            dy = A2 - A1;

            for ti=1:nlose
                logfif0 = log(mvnpdf(dy, mu0, sigma{ti})/mvnpdf(dy, mu0, sigma0));
                Wn(ti) = Wn(ti) + logfif0;
            end
            Wn = subplus(Wn);
            Wnvec(tao,:) = Wn;

            if sum(subplus(Wn-Tthreshold)) > 0 
                %if yes, update the average detection delay ADD
                if pathId==0
                    ADD = tao-1;
                else
                    ADD = (pathId/(pathId+1))*ADD + (tao-1)/(pathId+1);
                end

                lidx = find(Wn>Tthreshold);
                line_IDed(pathId+1,1:length(lidx)) = lidx;
                isolated_g0(pathId+1) = ~ismember(out,lidx);

                f_stop=1;

                break;
            end

            %Update the stopping variable till we actually stop
            tao=tao+1;
        end

        if mod(pathId,200) == 0
            pathId
            ADD
        end

    end
    
    if kk == 1 || (kk ~= 1 && ADD > ADDvec(kk-1))
        ADDvec(kk) = ADD;
        PFI_g0(kk) = sum(isolated_g0)/(nPath+1);
        kk = kk + 1;
    end

end
toc

save('L75out.mat', 'Avec', 'beta', 'ADDvec','PFI_g0','line_IDed');