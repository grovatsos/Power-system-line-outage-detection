clear all
close all

load('bus118info.mat')

n = size(B0,1);
m = length(B);

sigma = cell(m,1);

xsig = 0.03*eye(n);

%%-----------aaaaaa-------------------------
% sigma0 = B0*xsig.^2*B0';
% for i = 1:m
%     sigma{i} = B{i}*xsig.^2*transpose(B{i});
% end
%%-----------bbbbbb-------------------------
sigma0 = B0*2*xsig.^2*B0';
for i = 1:m
    sigma{i} = B{i}*2*xsig.^2*transpose(B{i});
end
%%------------------------------------------

x = xsig*randn(n,1);
mu0 = zeros(n,1);
A = log(n*1000);
nPath = 100;
Tthreshold = A*ones(1,m);
isolated_g0 = zeros(nPath+1,1);
line_IDed = zeros(nPath+1, 10);

for pathId=0:nPath
    
    tao=1;%stopping time
    f_stop=0;
    
    Wn=zeros(1,m);
    Wnvec = zeros(1e3,m);
    
    while f_stop == 0
        
%%------------aaaaaa------------------------
%         dx = xsig*randn(n,1);
%%------------bbbbbb------------------------
        xnew = xsig*randn(n,1);
        dx = xnew - x;
        x = xnew;
%%------------------------------------------
        
        dy = B0*dx;
        
        for ti=1:m
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
            
            f_stop=1;
            
            break;
        end
        
        %Update the stopping variable till we actually stop
        tao=tao+1;
    end
    
    pathId
    ADD
    
end