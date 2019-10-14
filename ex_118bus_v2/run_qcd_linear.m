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

out = 12;

beta = 500;
mu0 = zeros(nbus-1,1);
A = log(nlose*100);
nPath = 5000;
Tthreshold = A*ones(1,nlose);
isolated_g0 = zeros(nPath+1,1);
line_IDed = zeros(nPath+1, 10);

for pathId=0:nPath
    
    tao=1;%stopping time
    f_stop=0;
    
    Wn=zeros(1,nlose);
    Wnvec = zeros(1e3,nlose);
    
    while f_stop == 0

        x1 = Psig_*randn(nbus-1,1);
        x2 = Psig_*randn(nbus-1,1);
        dx = x2 - x1;

        dy = B0*dx;
        
        for ti=1:nlose
            logfif0 = log(mvnpdf(dy, mu0, sigma{ti})/mvnpdf(dy, mu0, sigma0));
            Wn(ti) = Wn(ti) + logfif0;
        end
        Wn = subplus(Wn);
        Wnvec(tao,:) = Wn;
        
        if sum(subplus(Wn-Tthreshold)) > 0 || tao > 3*beta
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
    
    pathId
    ADD
    
end