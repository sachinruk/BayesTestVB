function [mu,L,alpha_,beta,a,b]=vb3(y,x,alpha)
%performs variataional bayes and computes the best mu, Sigma etc.
%Input:
%y: response variable
%x: input variable: MUST have a column of ones if in form y=b0+b1*x
%alpha: asymmetricness: 0.5=symmetric

%first beta estimate
N=length(y); 
lga=0; lgb=0;
mu= x\y; 
%calculate approximate posterior mean
SigmaInv=(alpha*(1-alpha)/2)*(x'*x); 
%add a large variance to the diagonal- thanks to prior
D=size(SigmaInv,1); SigmaInv(1:(D+1):end)=SigmaInv(1:(D+1):end)+1/1000;
L=chol(SigmaInv)';
Lx=(L\x')';
gamma=-(1-2*alpha)/2*sum((y-x*mu));
delta=alpha*(1-alpha)/4*sum((y-x*mu).^2+sum(Lx.^2,2));
diff=0;

opts = optimset('MaxIter',1000,...
                'Diagnostics', 'off',...
                'GradObj', 'on',...
                'Display', 'off',...
                'TolFun',  1e-15,...
                'MaxFunEvals', 900);
            
for i=1:200
    [lga,lgb,l_one_sigma,l_one_sigma2]=q_sigma(lga,lgb,gamma,delta,N,opts);
    [alpha_,beta,l_one_w]=q_w(y,x,mu,Lx,l_one_sigma2,alpha);
    [mu, Lx,L,gamma,delta]=q_beta(y,x,l_one_sigma,l_one_sigma2,l_one_w,alpha);
    if i>1
        diff(i)=norm(params_old-[mu; lga; lgb; alpha_;beta]); 
        if diff(i)<1e-4, break; end
    end
    params_old=[mu; lga; lgb; alpha_;beta];
end
a=exp(lga); b=exp(lgb);

function [lg_a,lg_b,log_one_sigma,log_one_sigma2]=q_sigma(lg_a,lg_b,gamma,delta,N,opts)

par0=[lg_a; lg_b];            
[params, ~]=fminunc(@(params)F(params,gamma, delta,N),par0,opts);
lg_a=params(1); lg_b=params(2);                                                     
log_one_sigma=lg_a-lg_b; 
lg_a_1=log(exp(lg_a)+1); %log(a+1)
if lg_a>10, lg_a_1=lg_a; end
log_one_sigma2=log_one_sigma+lg_a_1-lg_b;

function [logF, derivative]=F(params,gamma, delta, N)
a=exp(params(1)); b=exp(params(2)); %need and b to be positive 
logF=(a-N)*(log(b)-psi(a))+(b-gamma)*a/b-delta*a*(a+1)/b^2 ...
                                                      -a*log(b)+gammaln(a);
logF=-logF;                                                  
partial_a=(N-a)*psi(1,a)-gamma/b-delta*(2*a+1)/b^2+1;                                                  
ln_partial_b=-N+gamma*a/b+2*delta*a*(a+1)/b^2;   
derivative=-[partial_a*a; ln_partial_b]; 


function [alpha_,beta,l_one_w]=q_w(y,x,mu,Lx,l_one_sigma2,alpha)
alpha_=(1-2*alpha)^2/(2*alpha*(1-alpha))+2;
ln_beta=l_one_sigma2+log(alpha)+log((1-alpha)/2)+...
                                          log((y-x*mu).^2+sum(Lx.^2,2));
beta=exp(ln_beta);
l_one_w=0.5*(log(alpha_)-ln_beta);


function [mu, Lx,L,gamma,delta]=q_beta(y,x,l_one_sigma,l_one_sigma2,l_one_w,alpha)
one_w_one_sigma2=exp(l_one_w+l_one_sigma2); one_w=exp(l_one_w);
%calculate approximate posterior mean
SigmaInv=(alpha*(1-alpha)/2)*bsxfun(@times,one_w_one_sigma2,x)'*x; 
%add a large variance to the diagonal- thanks to prior
D=size(SigmaInv,1); SigmaInv(1:(D+1):end)=SigmaInv(1:(D+1):end)+1/1000;
L=chol(SigmaInv)';
c1=exp(l_one_w+l_one_sigma2+log(alpha)+log((1-alpha)/2)); 
c2=(1-2*alpha)/2*exp(l_one_sigma);
mu=c1.*y-c2;
mu=sum(bsxfun(@times,mu,x))'; %sum over all N
mu=L'\(L\mu); 
gamma=-(1-2*alpha)/2*sum((y-x*mu));
Lx=(L\x')';
delta=alpha*(1-alpha)/4*sum(one_w.*((y-x*mu).^2+sum(Lx.^2,2)));
