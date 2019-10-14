clear all
close all

[bus line] = buslinedata(30);
line(29,2) = 22; %This is only for the 30 bus system to match diagram
nbus = size(bus,1);
nline = size(line,1);

A1 = [1:7,9:24];
A2 = [8,25:30];
interarea = ismember(line(:,1),A1)&ismember(line(:,2),A2) | ...
     ismember(line(:,1),A2)&ismember(line(:,2),A1);
interarea = find(interarea==1); %interarea contains the line #'s that are tie lines

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

J11 = J0(A1,A1);    J11(:,1) = [];  J11(1,:) = []; %bus 1 is assumed slack bus
J12 = J0(A1,A2);    J12(1,:) = [];
J21 = J0(A2,A1);   %J21 never used
J22 = J0(A2,A2);

S = get_line_flow(line, soln, ybus);
iflow = real(S(interarea));
interflow = zeros(nbus,1); %line flow assigned to "from" bus
for i = 1:length(iflow)
    interflow(line(interarea(i),1)) = interflow(line(interarea(i),1)) +...
        iflow(i);
end
interflow = interflow(A1); %internal area tie lines flows, assigned to from bus

J12_tilde = sum(J12,2); %row sum
J12_tilde = diag(J12_tilde);

J0(:,1) = [];   J0(1,:) = []; %remove slack bus

Psig = 0.03;

% theta0(1) =[];
% y = theta0;
% x = soln(2:end,4) - soln(2:end,6);
% 
% busk = bus;
% Pk = bus(2:end,4) + randn(nbus-1,1)*Psig;
% busk(2:end,4) = Pk;
% [soln ybus J] = get_pf(busk,line);
% S = get_line_flow(line, soln, ybus);
% iflownew = real(S(interarea));
% interflownew = zeros(nbus,1);
% for i = 1:length(iflownew)
%     interflownew(line(interarea(i),1)) = interflownew(line(interarea(i),1)) +...
%         iflownew(i);
% end
% interflownew = interflownew(A1);
% xnew = soln(2:end,4) - soln(2:end,6);
% ynew = soln(2:end,3);
% dx = xnew-x;
% dy = (ynew-y)*pi/180;
% dflow = interflownew(2:end)-interflow(2:end);
% x = xnew;
% y = ynew;
% interflow = interflownew;
% 
% busk = bus;
% Pk = bus(2:end,4) + randn(nbus-1,1)*Psig;
% busk(2:end,4) = Pk;
% [soln ybus J] = get_pf(busk,line);
% S = get_line_flow(line, soln, ybus);
% iflownew = real(S(interarea));
% interflownew = zeros(nbus,1);
% for i = 1:length(iflownew)
%     interflownew(line(interarea(i),1)) = interflownew(line(interarea(i),1)) +...
%         iflownew(i);
% end
% interflownew = interflownew(A1);
% xnew = soln(2:end,4) - soln(2:end,6);
% ynew = soln(2:end,3);
% dx = xnew-x;
% dy = (ynew-y)*pi/180;
% dflow = interflownew(2:end)-interflow(2:end);
% x = xnew;
% y = ynew;
% interflow = interflownew;
% 
% % % [-J12_tilde*dy(A1(2:end)-1)+J12*dy(A2-1) dflow]
% % % [dx(A1(2:end)-1)-dflow (J11+J12_tilde)*dy(A1(2:end)-1)]
% % % [dx(A1(2:end)-1) J11*dy(A1(2:end)-1)+J12*dy(A2-1) ...
% % %     J11*dy(A1(2:end)-1)+dflow+J12_tilde*dy(A1(2:end)-1)]
% % [dy(A1(2:end)-1) (J11+J12_tilde)\dx(A1(2:end)-1)-(J11+J12_tilde)\dflow]
% % figure
% % plot(1:length(ans), ans)
% % hold all; grid on
% % [dy J0\dx]
% % plot(1:length(ans), ans)

B0 = inv(J11+J12_tilde);
mu0 = zeros(length(A1)-1,1); %don't include slack
sigma0 = B0*(2*Psig^2*eye(length(A1)-1))*B0';

loseline = 1:nline;
temp = find(ismember(line(:,1),A2)|ismember(line(:,2),A2));
loseline(unique([13 16 temp'])) = [];
% loseline(unique([13 16 30 temp'])) = [];
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
    fij = fij(A1);
    fij(1) = [];
    dJ = -1/xij*(fij*fij');
    Ji = J11 + dJ + J12_tilde;
%     if rank(Ji) < length(Ji)
%         i
%     end
    Bi = inv(Ji);

    B{i} = Bi;
    sigma{i} = Ji\(2*Psig^2*eye(length(A1)-1))/(Ji');
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';

%     
% %     mu{i} = J0\((theta0(from)-theta0(to))*pi/180*(1/xij)*fij);
% %     
%     klmean(i) = 1/2*(mu{i}'*inv(sigma0)*mu{i});
    klcov(i) = 1/2 *(trace(sigma0\sigma{i}) - ...
        (length(A1)-1-1) + log(det(sigma0/sigma{i})));
end