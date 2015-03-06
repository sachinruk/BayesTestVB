clear all
% close all
clc

sigma=3;
alpha=0.25;

N=50000;
%richards simulations
u=rand(N,1); e=zeros(size(u));
e(u<alpha)=sigma/(1-alpha).*log(u(u<alpha)/alpha);
e(u>alpha)=sigma/alpha.*log((1-alpha)./(1-u(u>alpha)));

%my simulations
u=rand(N,1); w=-log(1-u);
z=randn(N,1);
e2=((1-2*alpha)/(alpha*(1-alpha)))*sigma*w+...
                        sqrt(2*w/(alpha*(1-alpha)))*sigma.*z;

%compare distributions                    
figure; subplot(121); hist(e,100); subplot(122); hist(e2,100)
figure; qqplot(e,e2)
disp([mean(e) mean(e2)]); 
disp([std(e) std(e2)]);
disp([skewness(e) skewness(e2)]);
disp([kurtosis(e) kurtosis(e2)]);