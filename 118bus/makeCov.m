function [sigma0 sigma] = makeCov(B0,B,xsig)

m = length(B);
sigma = cell(m,1);

sigma0 = B0*xsig^2*transpose(B0);
aa = triu(sigma0);
sigma0 = aa + triu(aa,1)';

for i = 1:m
    sigma{i} = B{i}*xsig^2*transpose(B{i});
    aa = triu(sigma{i});
    sigma{i} = aa + triu(aa,1)';
end

end