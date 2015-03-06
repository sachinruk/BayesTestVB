clear all
close all
% clc

beta=[13; 27]; d=length(beta);
N=100;
X=randn(N,d);
y=X*beta+randn(N,1);
t=double(y>0);
[mu, SigmaInv,C,Ey,Ey2]=vb_classification(t,X);
l=lb_classification(X,t,mu,SigmaInv,C);
disp(mu); disp(l);