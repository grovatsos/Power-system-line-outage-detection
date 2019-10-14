function F = fun3bus(x)

global X12 X23 X13 P2 P3 Q2 Q3

% F = [1/X12*(x(1)-x(2)) + 1/X13*(x(1)-x(3)) + 0.7;
%      1/X12*(x(2)-x(1)) + 1/X23*(x(2)-x(3)) - 0.2;
%      1/X13*(x(3)-x(1)) + 1/X23*(x(3)-x(2)) - 0.5];

% % flat voltage
% F = [1/X12*(x(1)) + 1/X23*(x(1)-x(2)) - P2;
%      1/X13*(x(2)) + 1/X23*(x(2)-x(1)) - P3];

% full power flow
t2 = x(1); t3 = x(2);
V2 = x(3); V3 = x(4);
F = [V2/X12*sin(t2) + V2*V3/X23*sin(t2-t3) - P2;
    V3/X13*sin(t3) + V2*V3/X23*sin(t3-t2) - P3;
    V2^2*(1/X12+1/X23) - V2/X12*cos(t2) ...
            - V2*V3/X23*cos(t2-t3) - Q2;
    V3^2*(1/X13+1/X23) - V3/X13*cos(t3) ...
            - V2*V3/X23*cos(t2-t3) - Q3];

end