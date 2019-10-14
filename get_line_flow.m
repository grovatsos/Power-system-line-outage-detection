function S = get_line_flow(line, soln, ybus)

nline = size(line,1);

from = line(:,1);
to = line(:,2);

Vmag = soln(:,2);
Vang = soln(:,3)*pi/180;
V = Vmag.*exp(1i*Vang);

% From the Ybus, extract the line impedances
yline = zeros(nline, 1);
for i = 1:nline
    yline(i) = -ybus(from(i), to(i));
end
I = yline .* (V(from) - V(to)); %current
S = V(from) .* conj(I); %power
% h = real(S);

end