%Testing AR1 Model


close all
clear all
seed=1
% rng(seed)

numsamples=100;
P=zeros(1, numsamples+1);
P(1)=0;
P0=0;
DeltaP=zeros(1, numsamples-1);
DeltaP2=zeros(1, numsamples/2);

for i=1:numsamples
    P(i+1)= P0+randn*0.5;
    DeltaP(i)=P(i+1)-P(i);
end

for i=1:numel(DeltaP2)
    DeltaP2(i)=P(2*i)-P(2*i-1);
end

figure
plot(P)
figure
plot(DeltaP)

figure
plot(DeltaP2)