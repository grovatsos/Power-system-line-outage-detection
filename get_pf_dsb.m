%distributed slack bus computation

function [soln, ybus, J, Jfull] = get_pf_dsb(mpc,alpha)
%values
eps = 1e-6;
nbus = size(mpc.bus,1);             %number of buses
nline = size(mpc.branch,1);           %number of transmission lines
V = mpc.bus(:,8);                   %bus voltage
ang = mpc.bus(:,9)*pi/180;            %phase (radians)
pt = 0;                         %real power losses in the system
Pgen = zeros(nbus,1);
Qgen = Pgen;
Pgen(mpc.gen(:,1),1)=mpc.gen(:,2);
Qgen(mpc.gen(:,1),1)=mpc.gen(:,3);
P = (Pgen-mpc.bus(:,3));               %bus power
Q = (Qgen-mpc.bus(:,4));               %bus reactive power
type = round(mpc.bus(:,2));        %type of bus
error = 1;                   
%error for while loop
count = 1;                      %count for while loop
v1 = find(type==3);             %angle ref bus index
v2 = zeros(nbus,1);             %pointing vector for type 2, PV buses
v3 = zeros(nbus,1);             %pointing vector for type 1, PQ buses
v4 = zeros(nbus,1);             %pointing vector for non-angle ref PV buses       
load = mpc.bus(:,3:4);     %initial value for load data

%assign pointing vector indices
a = find(type==2|type==3);
t2 = length(a);                 %number of PV buses
v2(1:t2) = a;
a = find(type==1);
t3 = length(a);                 %number of PQ buses
v3(1:t3) = a;
a = find(type==2);
t4 = length(a);                 %number of PV buses, sans the ref angle bus
v4(1:t4) = a;

%set initial voltage guesses
%V(v3(1:t3)) = ones(t3,1); inititial guess set above to nominal solution
v2o = v2;
t2o = t2;
v3o = v3;
t3o = t3;

%check and set generator limits
if size(mpc.gen,2)>11
    gmax = mpc.gen(:,4);
    gmin = mpc.gen(:,5);
else
    gmax = 999*ones(nbus,1);
    gmin = -999*ones(nbus,1);
end

%values for the ybus
link = mpc.branch(:,1:2);
phase = mpc.branch(:,10);
ratio = mpc.branch(:,9)+eq(mpc.branch(:,9),0);
tap = ratio.*exp(1j*phase*pi/180);
y = 1./(mpc.branch(:,3)+1i*mpc.branch(:,4));
ybus = zeros(nbus);

%compute the ybus
for k=1:nline
    w = [(y(k)+1i*mpc.branch(k,5)/2)/(tap(k)*tap(k)) -y(k)/conj(tap(k));
         -y(k)/tap(k) (y(k)+1i*mpc.branch(k,5)/2)];
    
    ybus(link(k,:),link(k,:)) = ybus(link(k,:),link(k,:))+w;
end

ybus = ybus+diag(mpc.bus(:,5)+1i*mpc.bus(:,6));

while((count < 20)&&(error > 10^-6))
    %allocate space for simulink
    S = complex(zeros(nbus,1));
    % 2 more rows/cols in J
    % J = complex(zeros(2*nbus));
    J = complex(zeros(t2+2*t3));

    x = zeros(2*nbus,1);
    flag = 0;

    %complex power
    for m=1:nbus
        for n=1:nbus
            S(m) = S(m)+conj(ybus(m,n))*V(m)*V(n)*exp(1i*(ang(m)-ang(n)));
        end
    end
    
    %check reactive limits
    for k=1:t2o
        u = gt(imag(S(v2o(k))),gmax(mpc.gen(:,1)==v2o(k))+10^-6);
        l = lt(imag(S(v2o(k))),gmin(mpc.gen(:,1)==v2o(k))-10^-6);
        if (u||l)&&eq(t2o,t2)
            Q(v2o(k)) = u*gmax(v2o(k))+l*gmin(v2o(k));
            S(v2o(k)) = P(v2o(k))+1i*Q(v2o(k));
            type(v2o(k)) = 1;
            flag = 1;
        end
    end
    
    %change index vectors if there is a violation
    if flag==1
        v2 = zeros(nbus,1);      
        v3 = zeros(nbus,1);            
        a = find(type==2|type==3);
        t2 = length(a);              
        v2(1:t2) = a;
        a = find(type==1);
        t3 = length(a);              
        v3(1:t3) = a;
        a = find(type==2);
        t4 = length(a);
        v4(1:t4) = a;
    end
    
    %mismatch equations
    
    fp = real(S)-P+pt*alpha;
    fq = imag(S)-Q;

    %compute the jacobian
    %|H N M||theta|
    %|K L R||  V  |
    %       |  pt | %this is the real power mismatch
    H = complex(zeros(nbus));
    N = complex(zeros(nbus));

    %dP/dtheta and dQ/dtheta
    %dP/dV and dQ/dV
    for m=1:nbus
        for n=1:nbus
            if(m == n)
                H(m,m) = H(m,m);
                N(m,m) = N(m,m)+2*conj(ybus(m,m))*V(m);
            else
                H(m,m) = H(m,m)+1i*conj(ybus(m,n))*V(m)*V(n)*exp(1i*(ang(m)-ang(n)));
                H(m,n) = -1i*conj(ybus(m,n))*V(m)*V(n)*exp(1i*(ang(m)-ang(n)));

                N(m,m) = N(m,m)+conj(ybus(m,n))*V(n)*exp(1i*(ang(m)-ang(n)));
                N(m,n) = conj(ybus(m,n))*V(m)*exp(1i*(ang(m)-ang(n)));
            end        
        end
    end

    K = imag(H);
    H = real(H);
    L = imag(N);
    N = real(N);
    M = alpha;
    R = zeros(t3,1);
    
    Jfull = [[H N M];[K L zeros(nbus,1)]];
    %make jacobian
    J(1:t2+2*t3,1:t2+2*t3) =...
        [H([v2(1:t2); v3(1:t3)],[v4(1:t4); v3(1:t3)]) ...
         N([v2(1:t2); v3(1:t3)],v3(1:t3)) ...
         M;
         K(v3(1:t3),[v4(1:t4); v3(1:t3)]) ...
         L(v3(1:t3),v3(1:t3)) ...
         R];

    %compute new values
    x(1:t2+2*t3) = [ang([v4(1:t4); v3(1:t3)]); V(v3(1:t3)); pt(1)]...
                    -real(J(1:(t2+2*t3),1:(t2+2*t3)))\...
                    [fp([v2(1:t2); v3(1:t3)]); fq(v3(1:t3))];
    
    %compute error
    error = norm(x(1:t2+2*t3)-[ang([v4(1:t4); v3(1:t3)]); V(v3(1:t3)); pt]);
    
    %assign new values
    ang([v4(1:t4); v3(1:t3)]) = x(1:t4+t3);
    V(v3(1:t3)) = x(t4+t3+1:t4+2*t3);
    pt(1) = x(t4+2*t3+1,1);
    count = count+1;
end   
S(abs(real(S))<eps) = 0;
%output
soln = zeros(nbus,13);
soln(:,1:3) = [mpc.bus(:,1) V ang*180/pi];
soln([v2o(1:t2o)],4:5) =  [-pt*alpha(v2o(1:t2o))+Pgen(v2o(1:t2o))...%-P(v2o(1:t2o))+real(S(v2o(1:t2o)))...
                             imag(S([v2o(1:t2o)]))];
% soln(v3o(1:t3o),6:7) =  -[real(S(v3o(1:t3o))) imag(S(v3o(1:t3o)))];
soln(:,6:7) =  [soln(:,4)-real(S) -imag(S).*eq(type,1)];

%soln(:,10:12) = [type gmax gmin];
soln(:,13) = pt*alpha;