function [J0 B0 B loseline] = makeB(mpc)

mpopt = mpoption('PF_ALG', 1, 'VERBOSE', 0, 'OUT_ALL', 0);
% [baseMVA, bus, gen, line, success, et, ybus, Yf, Yt, J] = 

results = runpf(mpc, mpopt);
baseMVA = results.baseMVA;
bus = results.bus;
line = results.branch;
[Ybus, ~, ~] = makeYbus(baseMVA, bus, line);

nbus = size(bus,1);
nline = size(line,1);
type = round(bus(:,2));

Ybus(type==3,:) = [];
Ybus(:,type==3) = [];

%% Aggregate multiple lines between the same 2 buses
[u,I,~] = unique(line(:,1:2), 'rows');
ixDupRows = setdiff(1:size(line(:,1:2),1), I);

%% Build base-case Jacobian matrix
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
J0(:,type==3) = [];
J0(type==3,:) = [];

%% Compute B matrices for base case and outage cases
B0 = J0\eye(length(J0)); %inv J

loseline = zeros(nline,1);
B = cell(nline,1);
cnt = 0;
for i = 1:nline
    
    from = line(i,1);
    to = line(i,2);
    xij = line(i,4);
    fij = zeros(nbus,1);
    fij(from) = 1;
    fij(to) = -1;
    fij(type==3) = [];
    dJ = -1/xij*(fij*fij');
    Ji = J0 + dJ;
    
    if rank(Ji) == size(Ji,1)
        cnt = cnt + 1;
        loseline(cnt) = i;
        B{cnt} = Ji\eye(length(Ji));
    end

end

B = B(1:cnt);
loseline = loseline(1:cnt);

end