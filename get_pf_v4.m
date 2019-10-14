% Nonlinear Newton-Rhapson Power Flow version 2
% Power System Toolbox test code
%%%%%%%%%%%%%%%%%%%%%%%%
% bus_sol = loadflow(bus,line,1e-9,10, 1.0,'n',1); %Run this line of code for obtaining the solution in PST 
%%%%%%%%%%%%%%%%%%%%%%%%
% Bus data format
    % col1 number
    % col2 voltage magnitude(pu)
    % col3 voltage angle(degree)
    % col4 p_gen(pu)
    % col5 q_gen(pu),
    % col6 p_load(pu)
    % col7 q_load(pu)
    % col8 G shunt(pu)
    % col9 B shunt(pu)
    % col10 bus_type
    %       bus_type - 1, swing bus
    %                - 2, generator bus (PV bus)
    %                - 3, load bus (PQ bus)
    % col11 q_gen_max(pu)
    % col12 q_gen_min(pu)
% Line data format
    % col1 from bus
    % col2 to bus
    % col3 resistance(pu)
    % col4 reactance(pu)
    % col5 line charging(pu)
    % col6 tap ratio
    % col7 phase shift(deg)
    
function [soln ybus dVdQ H K N L] = get_pf_v4(bus,line,tol,iter)

% default values when not defined
if ~exist('tol','var') 
    tol = 10^-6;
end

if ~exist('iter','var')
    iter = 25;
end
 
% values
nbus      = size(bus,1);                                % number of buses
typ       = bus(:,10);                                  % bus types
flag      = false;                                      % generator limit error
V         = bus(:,2);                                   % voltage magnitude
V(typ==3) = 1;                                          % reset inital guesses for PQ
ang       = [bus(typ==1,3)*pi/180; zeros(nbus-1,1)];    % angle converted (radians)
P         = bus(:,4)-bus(:,6);                          % bus power
Q         = bus(:,5)-bus(:,7);                          % bus reactive power
error     = 1;                                          % mismatch error
count     = 0;                                          % number of iterations
soln      = bus;                                        % initialize solution to bus

%check and set generator limits
if (size(bus,2)>11)
    gmax = bus(:,11);
    gmin = bus(:,12);
else
    gmax = 999*ones(nbus,1);
    gmin = -999*ones(nbus,1);
    soln(:,11:12) = [gmax,gmin];
end

% values for the ybus
link  = line(:,1:2);
phase = line(:,7);
ratio = line(:,6)+eq(line(:,6),0);
tap   = ratio.*exp(1j*phase*pi/180);
y     = 1./(line(:,3)+1i*line(:,4));
ybus  = zeros(nbus);

% compute the ybus
for k=1:size(line,1)
    w = [(y(k)+1i*line(k,5)/2)/(tap(k)*tap(k)) -y(k)/conj(tap(k));
         -y(k)/tap(k) (y(k)+1i*line(k,5)/2)];
    ybus(link(k,:),link(k,:)) = ybus(link(k,:),link(k,:))+w;
end
ybus = ybus+diag(bus(:,8)+1i*bus(:,9));

while (count<iter)&&((error>tol)||flag)
    % complex power and partial derivatives
    % |P| = |H N||theta|
    % |Q| = |K L||  V  |
    S  = conj(ybus*(V.*exp(1i*ang))).*V.*exp(1i*ang);
    fp = real(S)-P;
    fq = imag(S)-Q;
    
    H  = -1i*diag(V.*exp(1i*ang))*conj(ybus*diag(V.*exp(1i*ang)))+... %n~=k case, see Table 6.5
         diag(1i*conj(ybus*(V.*exp(1i*ang))).*V.*exp(1i*ang)); %see page 316 of power system book. this line is for n=k case, hence off diag elements are 0
    N  = diag(V.*exp(1i*ang))*conj(ybus*diag(exp(1i*ang)))+...
         diag(conj(ybus*(V.*exp(1i*ang))).*exp(1i*ang));
    
    % formulate Jacobian
    t = find(typ~=1);
    v = find(typ==3);
    K = imag(H(v,t));
    H = real(H(t,t));
    L = imag(N(v,v));
    N = real(N(t,v));
    J = [H N; K L];
    
    % solution and error check 
    x      = [ang(t); V(v)]-J\[fp(t); fq(v)];
    error  = norm(x-[ang(t); V(v)]);
    ang(t) = x(1:length(t));
    V(v)   = x(length(t)+1:end);
    S      = conj(ybus*(V.*exp(1i*ang))).*V.*exp(1i*ang);
    count  = count+1;
    
    % throw a generator limit error
    if error<tol
        u = imag(S(typ==2))>gmax(typ==2);
        l = imag(S(typ==2))<gmin(typ==2);
        if sum(u|l)>0
            k = find(typ==2);
            x = gmax(k).*u+gmin(k).*l;
            typ(k(u|l)) = 3;
            P(k(u|l))   = real(S(k(u|l)));
            Q(k(u|l))   = x(u|l);
            flag = true;
            disp('generator limit violation')
        elseif flag
            % reset bus types
            typ = bus(:,10);
            flag = false;
        end
    end
    
    if count==iter
        disp('solution did not converge')
    end
end

% outputs
% disp('number of iterations');
% disp(count);
soln(:,[2 3])  = [V ang*180/pi];
soln(typ~=3,5) = imag(S(typ~=3));
dVdQ           = inv(L-K*(H\N));