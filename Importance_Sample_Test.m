% code for testing importance sampling concept



clear all
close all
clc

mu0=0;
sigma0=1;
mu1=4;
sigma1=3;
% seed=5;
% rng(seed);


sample_size=1000000;
threshold=4;

x=zeros(1,sample_size); 
x_imp=zeros(1,sample_size); %use importance sampling
count=0;
count_imp=0;
temp=0;

p = normcdf(threshold,mu0,sigma0);
predicted=(1-p)*sample_size %theoretical number of samples we should get above threshold


for i=1:sample_size
    x(i)=normrnd(mu0,sigma0);
    x_imp(i)=normrnd(mu1,sigma1);
    
     if x(i)>=threshold
            count=count+1;
     end
     
     if x_imp(i)>=threshold
            temp=temp+normpdf(x_imp(i), mu0, sigma0)/normpdf(x_imp(i), mu1, sigma1);
            count_imp=count_imp+1;
     end
    
     if mod(i,10000)==0
         disp(['i = ', num2str(i)]);
     end
         
end
% actual_imp=(temp/count_imp)*sample_size;


disp(['predicted = ', num2str(predicted)]);
disp(['actual = ', num2str(count)]);
disp(['actual by importance sampling = ', num2str(temp)]);


% 
%  actual=count/sample_size
%  actual_imp=count_imp/sample_size