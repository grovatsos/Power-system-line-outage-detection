function [soln ybus J] = get_pf(bus,line)

%values
nbus = size(bus,1);             %number of buses
nline = size(line,1);           %number of transmission lines
V = bus(:,2);                   %bus voltage
ang = zeros(nbus,1);            %phase (radians)
P = bus(:,4)-bus(:,6);          %bus power
Q = bus(:,5)-bus(:,7);          %bus reactive power
type = round(bus(:,10));        %type of bus
error = 1;                      %error for while loop
count = 1;                      %count for while loop
v1 = find(type==1);             %slack bus index
v2 = zeros(nbus,1);             %pointing vector for type 2
v3 = zeros(nbus,1);             %pointing vector for type 3
load = bus(eq(type,3),6:7);     %initial value for load data

%assign pointing vector indices
a = find(type==2);
t2 = length(a);                 %number of PV buses
v2(1:t2) = a;
a = find(type==3);
t3 = length(a);                 %number of PQ buses
v3(1:t3) = a;

%set initial voltage guesses
V(v3(1:t3)) = ones(t3,1);
v2o = v2;
t2o = t2;
v3o = v3;
t3o = t3;

%check and set generator limits
if size(bus,2)>11
    gmax = bus(:,11);
    gmin = bus(:,12);
else
    gmax = 999*ones(nbus,1);
    gmin = -999*ones(nbus,1);
end

%values for the ybus
link = line(:,1:2);
phase = line(:,7);
ratio = line(:,6)+eq(line(:,6),0);
tap = ratio.*exp(1j*phase*pi/180);
y = 1./(line(:,3)+1i*line(:,4));
ybus = zeros(nbus);

%compute the ybus
for k=1:nline
    w = [(y(k)+1i*line(k,5)/2)/(tap(k)*tap(k)) -y(k)/conj(tap(k));
         -y(k)/tap(k) (y(k)+1i*line(k,5)/2)];
    
    ybus(link(k,:),link(k,:)) = ybus(link(k,:),link(k,:))+w;
end

ybus = ybus+diag(bus(:,8)+1i*bus(:,9));

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
    
    %check limits
    for k=1:t2o
        u = gt(imag(S(v2o(k))),gmax(v2o(k))+10^-6);
        l = lt(imag(S(v2o(k))),gmin(v2o(k))-10^-6);
        if (u||l)&&eq(t2o,t2)
            Q(v2o(k)) = u*gmax(v2o(k))+l*gmin(v2o(k));
            S(v2o(k)) = P(v2o(k))+1i*Q(v2o(k));
            type(v2o(k)) = 3;
            flag = 1;
        end
    end
    
    %change index vectors if there is a violation
    if flag==1
        v2 = zeros(nbus,1);      
        v3 = zeros(nbus,1);            
        a = find(type==2);
        t2 = length(a);              
        v2(1:t2) = a;
        a = find(type==3);
        t3 = length(a);              
        v3(1:t3) = a;
    end
    
    %mismatch equations
    fp = real(S)-P;
    fq = imag(S)-Q;

    %compute the jacobian
    %|H N||theta|
    %|K L||  V  |
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

    %make jacobian
    J(1:t2+2*t3,1:t2+2*t3) =...
        [H([v2(1:t2); v3(1:t3)],[v2(1:t2); v3(1:t3)])...
         N([v2(1:t2); v3(1:t3)],v3(1:t3));
         K(v3(1:t3),[v2(1:t2); v3(1:t3)])...
         L(v3(1:t3),v3(1:t3))];

    %compute new values
    x(1:t2+2*t3) = [ang([v2(1:t2); v3(1:t3)]); V(v3(1:t3))]...
                    -real(J(1:(t2+2*t3),1:(t2+2*t3)))\...
                    [fp([v2(1:t2); v3(1:t3)]); fq(v3(1:t3))];
    
    %compute error
    error = norm(x(1:t2+2*t3)-[ang([v2(1:t2); v3(1:t3)]); V(v3(1:t3))]);
    
    %assign new values
    ang([v2(1:t2); v3(1:t3)]) = x(1:t2+t3);
    V(v3(1:t3)) = x(t2+t3+1:t2+2*t3);
    
    count = count+1;
end   

%output
soln = zeros(nbus,12);
soln(:,1:3) = [bus(:,1) V ang*180/pi];
soln([v1; v2o(1:t2o)],4:5) =  [real(S([v1; v2o(1:t2o)]))...
                             imag(S([v1; v2o(1:t2o)]))];
soln(v3o(1:t3o),6:7) =  -[real(S(v3o(1:t3o))) imag(S(v3o(1:t3o)))];
soln(:,10:12) = [type gmax gmin];