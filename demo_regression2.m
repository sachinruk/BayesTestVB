clear all
close all
clc

sigmas=[0.1 0.1 2 2];
betas=[-2 13; -5 65; 2 5; 10 2];
alphas=0.04:0.02:0.96;
N=10000;
for k=1:4
    sigma= sigmas(k); 
    beta_true=betas(k,:);
    rng(1);
    for i=1:length(alphas)
        alpha=alphas(i);

        
        %simulate the errors
        u=rand(N,1); w=-log(1-u);
        z=randn(N,1);
        e=((1-2*alpha)/(alpha*(1-alpha)))*sigma*w+...
                                sqrt(2*w/(alpha*(1-alpha)))*sigma.*z;

        %simulate observations
        x=[ones(N,1) 10*randn(N,1)+5]; y=x*beta_true'+e;
        [mu,L,alpha_,beta,a,b]=vb3(y,x,alpha);
        lb(i)=lowerBound(y,mu,L,alpha_,beta,a,b,alpha);

        % [mu2,Sigma,alpha_2,beta2,a2,b2]=vb2(y,x,alpha);
        % F(y,mu2,Sigma,alpha_2,beta2,a2,b2,alpha)

        %metropolis
        step=0.01;
        [beta, AR]=metropolis(y,x,alpha, step);
        % figure; plot(beta(:,1),beta(:,2),'.');
        Sigma=cov(beta);
        llik(i)=logLikIS(y,x,alpha,beta_true',Sigma);
        if mod(i,10)==0, disp(i); end;
    end
    % likMetropolis(y,x,beta(:,1000:end),alpha,mu,Sigma)

    figure; plot(alphas,llik); hold on; plot(alphas,lb,'r')
    % figure; plot(alphas,llik-lb); 
end

% opt.logpdf=@(beta)lg_pby(beta,y,x,alpha);
% opt.V=inv(L*L');
% opt.Dims=2;
% opt.dsp=0;
% opt.Mmax=1e5;
% opt.div=max(20,ceil(opt.Mmax/1e6));
% opt.hgrd=1e2;
% opt.filename='test';
% beta2=AM([0;0],opt);