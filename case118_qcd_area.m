clear all
close all

%% Initialize case
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

% Internal and exteral areas defined
A1 = [1:45, 113:115, 117];
A2 = [46:112, 116, 118];

% Find interarea lines
interarea = ismember(line(:,1),A1)&ismember(line(:,2),A2) | ...
     ismember(line(:,1),A2)&ismember(line(:,2),A1);
interarea = find(interarea==1);

% Variance of power injection data
% Psig = diag(0.03*ones(nbus-1,1));
% Psig = diag(0.001*bus(type~=1,4)) + ...
%         diag(0.001*bus(type~=1,6)) + ...
%         diag(0.003*ones(nbus-1,1));
Psig = diag(0.03*ones(nbus,1));

% Aggregate multiple lines between the same 2 buses
[u,I,J] = unique(line(:,1:2), 'rows');
hasDuplicates = size(u,1) < size(line(:,1:2),1);
ixDupRows = setdiff(1:size(line(:,1:2),1), I);
dupRowValues = line(ixDupRows,1:2);

% Build base-case Jacobian matrix
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

%% Area 1 slack - break Jacobian into 4 portions based on internal 
%  and external systems
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

y = (bus(type~=3,VA)-results.bus(type==3,VA))*pi/180;
x = -bus(:,3)/baseMVA;
x(type~=1) = x(type~=1) + results.gen(:,2)/baseMVA;

%% Generate the first two measurements under base-case system
busk = bus;
Pk = bus(:,3) + Psig*baseMVA*randn(nbus,1);
busk(:,3) = Pk;
mpc.bus = busk;
results = runpf(mpc, mpopt);
iflownew = results.branch(interarea,PF)/baseMVA;
interflownew = zeros(nbus,1);
for i = 1:length(iflownew)
    interflownew(line(interarea(i),1)) = interflownew(line(interarea(i),1)) +...
        iflownew(i);
end
% interflownew = interflownew(A1_);
interflownew = interflownew(A1);  % Area 1 slack
interflownew(type==3) = [];    % Area 1 slack
xnew = -results.bus(:,3)/baseMVA;
xnew(type~=1) = xnew(type~=1) + results.gen(:,2)/baseMVA;
dx = xnew(type~=3) - x(type~=3);
ynew = (results.bus(type~=3,VA)-results.bus(type==3,VA))*pi/180;
dy = ynew - y;
dflow = interflownew - interflow;
x = xnew;
y = ynew;
interflow = interflownew;

busk = bus;
Pk = bus(:,3) + Psig*baseMVA*randn(nbus,1);
busk(:,3) = Pk;
mpc.bus = busk;
results = runpf(mpc, mpopt);
iflownew = results.branch(interarea,PF)/baseMVA;
interflownew = zeros(nbus,1);
for i = 1:length(iflownew)
    interflownew(line(interarea(i),1)) = interflownew(line(interarea(i),1)) +...
        iflownew(i);
end
% interflownew = interflownew(A1_);
interflownew = interflownew(A1);  % Area 1 slack
interflownew(type==3) = [];    % Area 1 slack
xnew = -results.bus(:,3)/baseMVA;
xnew(type~=1) = xnew(type~=1) + results.gen(:,2)/baseMVA;
dx = xnew(type~=3) - x(type~=3);
ynew = (results.bus(type~=3,VA)-results.bus(type==3,VA))*pi/180;
dy = ynew - y;
dflow = interflownew - interflow;
x = xnew;
y = ynew;
interflow = interflownew;

% Check that the approximations used match nonlinear solutions
figure;
plot(1:length(A1_), dx(A1_), 1:length(A1_), ...
    J11*dy(A1_)+J12*dy(A2_), 1:length(A1_), ...
    (J11+J12_tilde)*dy(A1_)+dflow);
grid on; hold on;
% figure;
% dxx = -imag(Ybus)*dy;
% plot(1:length(A1_), dx(A1_), 1:length(A1_), dxx(A1_));
% hold on; grid on;
% figure;
% plot(1:length(A1_), dx(A1_)-dflow, 1:length(A1_), ...
%     (J11+J12_tilde)*dy(A1_));
% figure;
% plot(1:length(A1_), (J11+J12_tilde)\(dx(A1_)-dflow), 1:length(A1_), ...
%     dy(A1_));
% figure;
% plot(1:length(A1_), dy(A1_), 1:length(A1_), ...
%     (J11+J12_tilde)\dx(A1_)-(J11+J12_tilde)\dflow);
% hold all; grid on
% plot(1:length(dy), dy, 1:length(dy), J0\dx);
% figure;
% plot(1:length(A1_), J12*dy(A2_), 1:length(A1_), ...
%     dflow+J12_tilde*dy(A1_));

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

% Compute D(fi || f0)
for i = 1:numel(loseline)
    
    from = line(loseline(i),1);
    to = line(loseline(i),2);
    xij = line(loseline(i),4);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
    fij(type==3) = [];    % Area 1 slack
    fij = fij(A1_);
%     fij = fij(A1);
%     fij(type==3) = [];    % Area 1 slack
    dJ = -1/xij*(fij*fij');
    Ji = J11 + dJ + J12_tilde;
    Bi = inv(Ji);

    B{i} = Bi;
    sigma{i} = Ji\(2*Psig(1:length(A1_),1:length(A1_)).^2*...
        eye(length(A1_)))/(Ji');
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';
%     
% %     mu{i} = J0\((theta0(from)-theta0(to))*pi/180*(1/xij)*fij);
% %     
%     klmean(i) = 1/2*(mu{i}'*inv(sigma0)*mu{i});
    klcov(i) = 1/2 *(trace(sigma0\sigma{i}) - ...
        length(A1_) + log(det(sigma0/sigma{i})));
end

% Compute D(fi || f0) - D(fi || fj)
klcovij = zeros(numel(loseline),numel(loseline));
for i = 1:numel(loseline)
    for j = 1:length(loseline)
        klcovij(j,i) = klcov(i) - ...
            1/2 *(trace(sigma{j}\sigma{i}) - ...
            length(A1_) + log(det(sigma{j}/sigma{i})));
    end
end

% Look at 4 line outages
out = 12; % min KL div
% out = 30;
% out = 14;
% out = 7;  % max KL div

% Consider these FAR
beta = [1 6 12 24 48 168];
Lbeta = length(loseline)*30*3600.*beta;
Avec = log(Lbeta);  % A threshold corresponding to FARs
Avec = Avec(1);
ADDvec = zeros(1,length(Avec));
PFI_g0 = zeros(1,length(Avec));

total_paths = 1;    % Set total # of sample paths
isolate_percent = zeros(length(loseline),1);
out_est_struct = cell(length(loseline),1);

Anum = 1;
while Anum <= length(Avec)
% for Anum = 1:length(Avec)
    
    A=Avec(Anum) %Threshold
    ADD = 0;
    isolated_g0 = zeros(total_paths+1,1);
    isolated_max = zeros(total_paths+1,1);
    line_IDed = zeros(total_paths+1, 10);

    line_post = line;
    line_post(loseline(out),:) = [];
    interarea_post = interarea;
    interarea_post(interarea_post > loseline(out)) = ...
        interarea_post(interarea_post > loseline(out)) - 1;

    %----------------- For testing purposes ------------------------
% %     [u,I,J] = unique(line_post(:,1:2), 'rows');
% %     hasDuplicates = size(u,1) < size(line_post(:,1:2),1);
% %     ixDupRows = setdiff(1:size(line_post(:,1:2),1), I);
% %     dupRowValues = line_post(ixDupRows,1:2);
% %     
% %     J0post = zeros(nbus,nbus);
% %     for i = 1:nline-1
% %         from = line_post(i,1);
% %         to = line_post(i,2);
% %         Xl = line_post(i,4);
% %         J0post(from,to) = -1/Xl;
% %         J0post(to,from) = -1/Xl;
% %         J0post(from,from) = J0post(from,from) + 1/Xl;
% %         J0post(to,to) = J0post(to,to) + 1/Xl;
% %     end
% %     for i = ixDupRows
% %         from = line_post(i,1);
% %         to = line_post(i,2);
% %         Xl = line_post(i,4);
% %         J0post(from,to) = J0post(from,to) - 1/Xl;
% %         J0post(to,from) = J0post(to,from) - 1/Xl;
% %     end
% %     J11post = J0post(A1,A1);
% %     J11post(A1==find(type==3),:) = [];
% %     J11post(:, A1==find(type==3)) = [];
% %     J12post = J0post(A1,A2);
% %     J12post(A1==find(type==3),:) = [];
% %     J22post = J0post(A2,A2);
% %     
% %     busk = bus;
% %     Pk = bus(:,3) + Psig*baseMVA*randn(nbus,1);
% %     busk(:,3) = Pk;
% %     mpc.bus = busk;
% %     mpc.branch = line_post;
% %     results = runpf(mpc, mpopt);
% %     iflownew = results.branch(interarea_post,PF)/baseMVA;
% %     interflownew = zeros(nbus,1);
% %     for i = 1:length(iflownew)
% %         interflownew(line_post(interarea_post(i),1)) = ...
% %             interflownew(line_post(interarea_post(i),1)) + iflownew(i);
% %     end
% %     % interflownew = interflownew(A1_);
% %     interflownew = interflownew(A1);  % Area 1 slack
% %     interflownew(type==3) = [];    % Area 1 slack
% %     xnew = -results.bus(:,3)/baseMVA;
% %     xnew(type~=1) = xnew(type~=1) + results.gen(:,2)/baseMVA;
% %     dx = xnew(type~=3) - x(type~=3);
% %     ynew = (results.bus(type~=3,VA)-results.bus(type==3,VA))*pi/180;
% %     dy = ynew - y;
% %     dflow = interflownew - interflow;
% %     x = xnew;
% %     y = ynew;
% %     interflow = interflownew;
% %     
% %     busk = bus;
% %     Pk = bus(:,3) + Psig*baseMVA*randn(nbus,1);
% %     busk(:,3) = Pk;
% %     mpc.bus = busk;
% %     mpc.branch = line_post;
% %     results = runpf(mpc, mpopt);
% %     iflownew = results.branch(interarea_post,PF)/baseMVA;
% %     interflownew = zeros(nbus,1);
% %     for i = 1:length(iflownew)
% %         interflownew(line_post(interarea_post(i),1)) = ...
% %             interflownew(line_post(interarea_post(i),1)) + iflownew(i);
% %     end
% %     % interflownew = interflownew(A1_);
% %     interflownew = interflownew(A1);  % Area 1 slack
% %     interflownew(type==3) = [];    % Area 1 slack
% %     xnew = -results.bus(:,3)/baseMVA;
% %     xnew(type~=1) = xnew(type~=1) + results.gen(:,2)/baseMVA;
% %     dx = xnew(type~=3) - x(type~=3);
% %     ynew = (results.bus(type~=3,VA)-results.bus(type==3,VA))*pi/180;
% %     dy = ynew - y;
% %     dflow = interflownew - interflow;
% %     x = xnew;
% %     y = ynew;
% %     interflow = interflownew;
% %     
% %     [Ybuspost, Yf, Yt] = makeYbus(baseMVA, bus, line_post);
% %     Ybuspost(type==3,:) = [];
% %     Ybuspost(:,type==3) = [];
% %     figure;
% %     plot(1:length(A1_), dx(A1_), 1:length(A1_), ...
% %         J11post*dy(A1_)+J12post*dy(A2_), 1:length(A1_), ...
% %         (J11post+J12_tilde)*dy(A1_)+dflow);
% %     grid on; hold on;
% %     figure;
% %     dxx = -imag(Ybuspost)*dy;
% %     plot(1:length(A1_), dx(A1_), 1:length(A1_), dxx(A1_));
% %     hold on; grid on;
% %     figure;
% %     plot(1:length(A1_), dx(A1_)-dflow, 1:length(A1_), ...
% %         (J11post+J12_tilde)*dy(A1_));
% %     figure;
% %     plot(1:length(A1_), (J11post+J12_tilde)\(dx(A1_)-dflow), 1:length(A1_), ...
% %         dy(A1_));
% %     figure;
% %     plot(1:length(A1_), dy(A1_), 1:length(A1_), ...
% %         (J11+J12_tilde)\dx(A1_)-(J11+J12_tilde)\dflow);
% %     hold all; grid on
% %     plot(1:length(dy), dy, 1:length(dy), J0\dx);
% %     figure;
% %     plot(1:length(A1_), J12post*dy(A2_), 1:length(A1_), ...
% %         dflow+J12_tilde*dy(A1_));

    Tthreshold = A*ones(1,length(loseline));
    changetime = 1;
    
    %Run for 1000 sample paths
    pathIdx = 0;
    tic

    for pathId = 0:total_paths
        tao=1;%stopping time
        f_stop=0;
        Wn=zeros(1,numel(loseline));
        Wnvec = zeros(1e3,length(loseline));
        linek = line;
        interareak = interarea;
        bb = B0;

        while f_stop==0

            if tao == changetime
                linek = line_post;
                interareak = interarea_post;
                bb = B{out};
            end

            busk = bus;
            Pk = bus(:,3) + Psig*baseMVA*randn(nbus,1);
            busk(:,3) = Pk;
            mpc.branch = linek;
            mpc.bus = busk;
            results = runpf(mpc, mpopt);
            iflownew = results.branch(interareak,PF)/baseMVA;
            interflownew = zeros(nbus,1);
            for i = 1:length(iflownew)
                interflownew(linek(interareak(i),1)) = ...
                    interflownew(linek(interareak(i),1)) + iflownew(i);
            end
    %         interflownew = interflownew(A1_);
            interflownew(type==3) = [];        % Area 1 slack
            interflownew = interflownew(A1_);     % Area 1 slack
            xnew = -results.bus(:,3)/baseMVA;
            xnew(type~=1) = xnew(type~=1) + results.gen(:,2)/baseMVA;
            dx = xnew(type~=3) - x(type~=3);
            ynew = (results.bus(type~=3,VA)-results.bus(type==3,VA))*pi/180;
            dy_ = ynew - y;
            dflow = interflownew - interflow;
            x = xnew;
            y = ynew;
            interflow = interflownew;

    %         for ti=1:length(loseline)
    %            Wn(ti) = Wn(ti) + log(mvnpdf(dy, mu0, sigma{ti})/mvnpdf(dy, mu0, sigma0));
    %            if log(mvnpdf(dy, mu{ti}, sigma0)/mvnpdf(dy, mu0, sigma0)) > Wn(ti)
    %                Wn(ti) = log(mvnpdf(dy, mu{ti}, sigma0)/mvnpdf(dy, mu0, sigma0));
    %            end
    %         end
    %         Wn = subplus(Wn);

%             figure;
% %             plot(1:length(A1_), dy_(A1_), 1:length(A1_), ...
% %                     bb*dx(A1_)-bb*dflow);
%             plot(1:length(A1_), dx(A1_)-dflow, 1:length(A1_), ...
%                 bb\dy_(A1_));
%             grid on; hold on;

            for ti=1:length(loseline)
                dy = dy_(A1_) + B{ti}*dflow;
                mu0 = B{ti}*dflow;
                Wn(ti) = Wn(ti) + log(mvnpdf(dy, mu0, sigma{ti})/mvnpdf(dy, mu0, sigma0));
            end
            Wn = subplus(Wn);

            Wnvec(tao,:) = Wn;

            %Check if any one of them is above A
            if sum(subplus(Wn-Tthreshold))> 0
                lidx = find(Wn>Tthreshold);
                line_IDed(pathId+1,1:length(lidx)) = lidx;
                isolated_g0(pathId+1) = ~ismember(out,lidx);
%                 isolated_max(pathId+1) = find(Wn-Tthreshold == max(Wn-Tthreshold));
                
                if tao < changetime
                    if pathId==0 
                        PFA = 1;
                    else
                        PFA = (pathId/(pathId+1))*PFA + (1)/(pathId+1);
                    end
                else
                    if pathId==0 
                        PFA = 0;
                    else
                       PFA = (pathId/(pathId+1))*PFA + (0)/(pathId+1);
                    end

                    if pathIdx==0
                        ADD = tao-changetime;
                    else
                        ADD = (pathIdx/(pathIdx+1))*ADD + (tao-changetime)/(pathIdx+1);
                    end
                    pathIdx = pathIdx + 1;
                end

                f_stop=1;

                break;
            end

            %Update the stopping variable till we actually stop
            tao=tao+1;
        end
    end
    
    toc
    
    ADD
    
    if Anum == 1 || (Anum ~= 1 && ADD > ADDvec(Anum-1))
        ADDvec(Anum) = ADD;
    %     PFI_max(Anum) = sum(isolated_max~=out)/(total_paths+1);
        PFI_g0(Anum) = sum(isolated_g0)/(total_paths+1);
        Anum = Anum + 1;
    end
    
%     figure;
%     plot(Wnvec(1:tao,:));
%     hold on; grid on;
%     plot(Wnvec(1:tao,out), 'LineWidth',2);
    
    pause(1);
    
    % Generate 2 more measurements sets under the base-case
    % system - 'reset' the incremental differences for next
    % sample path
    busk = bus;
    Pk = bus(:,3) + Psig*baseMVA*randn(nbus,1);
    busk(:,3) = Pk;
    mpc.bus = busk;
    mpc.branch = line;
    results = runpf(mpc, mpopt);
    iflownew = results.branch(interarea,PF)/baseMVA;
    interflownew = zeros(nbus,1);
    for i = 1:length(iflownew)
        interflownew(line(interarea(i),1)) = interflownew(line(interarea(i),1)) +...
            iflownew(i);
    end
    % interflownew = interflownew(A1_);
    interflownew(type==3) = [];        % Area 1 slack
    interflownew = interflownew(A1_);     % Area 1 slack
    xnew = -results.bus(:,3)/baseMVA;
    xnew(type~=1) = xnew(type~=1) + results.gen(:,2)/baseMVA;
    dx = xnew(type~=3) - x(type~=3);
    ynew = (results.bus(type~=3,VA)-results.bus(type==3,VA))*pi/180;
    dy = ynew - y;
    dflow = interflownew - interflow;
    x = xnew;
    y = ynew;
    interflow = interflownew;

    busk = bus;
    Pk = bus(:,3) + Psig*baseMVA*randn(nbus,1);
    busk(:,3) = Pk;
    mpc.bus = busk;
    mpc.branch = line;
    results = runpf(mpc, mpopt);
    iflownew = results.branch(interarea,PF)/baseMVA;
    interflownew = zeros(nbus,1);
    for i = 1:length(iflownew)
        interflownew(line(interarea(i),1)) = interflownew(line(interarea(i),1)) +...
            iflownew(i);
    end
    % interflownew = interflownew(A1_);
    interflownew(type==3) = [];        % Area 1 slack
    interflownew = interflownew(A1_);     % Area 1 slack
    xnew = -results.bus(:,3)/baseMVA;
    xnew(type~=1) = xnew(type~=1) + results.gen(:,2)/baseMVA;
    dx = xnew(type~=3) - x(type~=3);
    ynew = (results.bus(type~=3,VA)-results.bus(type==3,VA))*pi/180;
    dy = ynew - y;
    dflow = interflownew - interflow;
    x = xnew;
    y = ynew;
    interflow = interflownew;

end

% save('L12out.mat', 'Avec', 'beta', 'ADDvec','PFI_g0','line_IDed');


% % ADDvec = zeros(length(loseline),1);
% % PFAvec = zeros(length(loseline),1);
% % total_paths = 19;
% % isolate_percent = zeros(length(loseline),1);
% % out_est_struct = cell(length(loseline),1);
% % for out_idx = 1:length(loseline)
% %     out_idx
% %     out = loseline(out_idx);
% %     line_post = line;
% %     line_post(out,:) = [];
% % 
% %     [soln ybus J] = get_pf(bus,line_post);
% %     V0_post = soln(1:end,2);
% %     theta0_post = soln(1:end,3);
% % 
% %     %nline;
% %     A = 1000; %Threshold
% %     ADD = 0;
% %     Tthreshold = A*ones(1,length(loseline));
% %     changetime = 1;
% %     out_est = zeros(total_paths+1,5);
% % 
% %     %Run for 1000 sample paths
% %     pathIdx = 0;
% %     tic
% %     isolate_cnt = 0;
% % 
% %     for pathId=0:total_paths
% %         tao=1;%stopping time
% %         f_stop=0;
% %         Wn=zeros(1,numel(loseline));
% %         Wnvec = zeros(1e3,length(loseline));
% %         linek = line;
% %         while f_stop==0
% % 
% %             %This statement will be removed and the X vector must come from
% %             %Christine's programs. May be a function call?
% %     %         dy=mvnrnd(mu, sigma0, 1);
% % 
% %             if tao >= changetime
% %                 linek = line_post;
% %             end
% % 
% %             busk = bus;
% %             Pk = bus(type~=1,4) + Psig*randn(nbus-1,1);
% %             busk(type~=1,4) = Pk;
% % 
% %             [soln ybus J] = get_pf(busk,linek);
% %             ynew = soln(type~=1,3);
% %             dy = (ynew-y)*pi/180;
% %             y = ynew;
% % 
% %             for ti=1:length(loseline)
% %                Wn(ti) = max(subplus(Wn(ti) + ...
% %                         log(mvnpdf(dy, mu0, sigma{ti})/mvnpdf(dy, mu0, sigma0))), ...
% %                         log(mvnpdf(dy,mu{ti},sigma0)/mvnpdf(dy, mu0, sigma0)));       
% %     %            Wn(ti) = Wn(ti) + log(mvnpdf(dy, mu0, sigma{ti})/mvnpdf(dy, mu0, sigma0));
% %     %            if log(mvnpdf(dy, mu{ti}, sigma0)/mvnpdf(dy, mu0, sigma0)) > Wn(ti)
% %     %                Wn(ti) = log(mvnpdf(dy, mu{ti}, sigma0)/mvnpdf(dy, mu0, sigma0));
% %     %            end
% %             end
% % 
% %     %         for ti=1:length(loseline)
% %     %             Wn(ti) = Wn(ti) + log(mvnpdf(dy, mu0, sigma{ti})/mvnpdf(dy, mu0, sigma0));
% %     %         end
% %     %         Wn = subplus(Wn);
% % 
% %             Wnvec(tao,:) = Wn;
% % 
% %             %Check if any one of them is above A
% %             if sum(subplus(Wn-Tthreshold))> 0
% %                 pathId;
% %                 tao;
% %                 aa = find(Wn>Tthreshold);
% %                 out_est(pathId+1,1:length(aa)) = aa;
% %                 if ismember(out_idx,aa)
% %                     isolate_cnt = isolate_cnt + 1;
% %                 end
% %                 if tao < changetime
% %                     if pathId==0 
% %                         PFA = 1;
% %                     else
% %                         PFA = (pathId/(pathId+1))*PFA + (1)/(pathId+1);
% %                     end
% %                 else
% %                     if pathId==0 
% %                         PFA = 0;
% %                     else
% %                        PFA = (pathId/(pathId+1))*PFA + (0)/(pathId+1);
% %                     end
% % 
% %                     if pathIdx==0
% %                         ADD = tao-changetime;
% %                     else
% %                         ADD = (pathIdx/(pathIdx+1))*ADD + (tao-changetime)/(pathIdx+1);
% %                     end
% %                     pathIdx = pathIdx + 1;
% %                 end
% % 
% %                 f_stop=1;
% % 
% %                 break;
% %             end
% % 
% %             %Update the stopping variable till we actually stop
% %             tao=tao+1;
% %         end
% %     end
% %     toc
% %     out_est
% %     out_est_struct{out_idx} = out_est;
% %     save('out_est_struct', 'out_est_struct');
% %     
% %     isolate_percent(out_idx) = isolate_cnt/(total_paths+1);
% %     ADDvec(out_idx) = ADD;
% %     PFAvec(out_idx) = PFA;
% % 
% % %     figure;
% % %     plot(Wnvec(1:tao,:));
% % %     grid on; hold on;
% % %     plot(Wnvec(1:tao,out_idx), 'LineWidth', 4);
% % end
    