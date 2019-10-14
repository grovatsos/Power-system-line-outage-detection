%Compute the KL divergence D(f1||f0) 
function [KL] = KL_compute(sigma1, sigma0, mu1, mu0)
    
    if nargin==2
        KL=1/2 *(trace(sigma0\sigma1)-(length(sigma0)) + log(det(sigma0/sigma1)));
    else
        KL=1/2 *((mu1-mu0)'*(sigma0\(mu1-mu0))+trace(sigma0\sigma1)-(length(sigma0)) + log(det(sigma0/sigma1)));
    end
end