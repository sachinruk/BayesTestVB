clear all
close all
% clc

a0=0.1;a1=0.1;b1=0.85;
n=2500;

t=zeros(1,5000);
rng(1); %set seed
for i=1:5
%     tstart=tic;
    [y, sig, VaR5, VaR1, VaR5HS, VaR1HS]=gengarchnHS2(n,a0,a1,b1);
    [BFucn1s2500(i),BFcc1n1s2500(i),BFind1n1s2500(i),BFcc4n1s2500(i),...
        BFind4n1s2500(i),BFvqn1s2500(i)]=runBFtests(y',VaR1',0.01);
%     BFdq2n1s2500(i) = 1/bfdq2mlik2(y,VaR1,0.01);
%     BFdq2n5s2500(i)=1/bfdq2mlik2(y,VaR5,0.05);
%     BFdq2n1p2500(i) = 1/bfdq2mlik2(y,VaR1HS,0.01);
%     BFdq2n5p2500(i)=1/bfdq2mlik2(y,VaR5HS,0.05);
  
%     t(i)=toc(tstart);
    disp(i);
end
%     [BFucn5s2500(i),BFcc1n5s2500(i),BFind1n5s2500(i),BFcc4n5s2500(i),BFind4n5s2500(i),BFvqn5s2500(i)]=runBFtests(y,VaR5,0.05);
%     [BFucn1p2500(i),BFcc1n1p2500(i),BFind1n1p2500(i),BFcc4n1p2500(i),BFind4n1p2500(i),BFvqn1p2500(i)]=runBFtests(y,VaR1HS,0.01);
%     [BFucn5p2500(i),BFcc1n5p2500(i),BFind1n5p2500(i),BFcc4n5p2500(i),BFind4n5p2500(i),BFvqn5p2500(i)]=runBFtests(y,VaR5HS,0.05);