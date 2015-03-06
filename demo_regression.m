clear all
close all
clc

sigma= 2.3; alpha=0.05;
beta=[1 2];

N=1000;
%simulate the errors
u=rand(N,1); w=-log(1-u);
z=randn(N,1);
e=((1-2*alpha)/(alpha*(1-alpha)))*sigma*w+...
                        sqrt(2*w/(alpha*(1-alpha)))*sigma.*z;

%simulate observations
x=[ones(N,1) 10*randn(N,1)+5]; y=x*beta'+e;
mu=vb(y,x,0.05);