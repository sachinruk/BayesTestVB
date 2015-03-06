function [mu,Sigma,alpha_,beta,a,b]=vb2(y,x,alpha)
%performs variataional bayes and computes the best mu, Sigma etc.
%Uses the SVD- bit slower than choleski
%Input:
%y: response variable
%x: input variable: MUST have a column of ones if in form y=b0+b1*x
%alpha: asymmetricness: 0.5=symmetric

%first beta estimate
N=length(y); 
l_one_w=log(ones(N,1)); a=1; b=1;
mu= x\y; 
% SigmaInv0=(alpha*(1-alpha)/2)*(x'*x); step=length(SigmaInv0);
% SigmaInv0(1:(step+1):end)=SigmaInv0(1:(step+1):end)+1/1000;
[U,D,~]=svd((bsxfun(@times,exp(l_one_w),x)'*x)); d=(alpha*(1-alpha)/2)*diag(D);
sigmaInv_d=d*a*(a+1)/b^2+1/1000; L=diagMult(U,1./sqrt(sigmaInv_d),'r');
traceTerm=x*L;
% L=chol(SigmaInv0);
gamma=-(1-2*alpha)/2*sum((y-x*mu));
% betaX=x*U;
delta=alpha*(1-alpha)/4*sum((y-x*mu).^2+sum(traceTerm.^2,2));
for i=1:200
%     [a,b,l_one_sigma,l_one_sigma2]=q_sigma(y,x,mu,l_one_w,...
%                                                     alpha,a,b,gamma,delta);
    [a,b,l_one_sigma,l_one_sigma2]=q_sigma(a,b,gamma,delta,N);
    [alpha_,beta,l_one_w]=q_w(y,x,mu,L,l_one_sigma2,alpha);
    [mu, L,gamma,delta]=q_beta(y,x,l_one_sigma,l_one_sigma2,l_one_w,alpha);
end
Sigma=L*L';

function [a,b,log_one_sigma,log_one_sigma2]=q_sigma(a,b,gamma,delta,N)
opts = optimset('MaxIter',1000,...
                'Diagnostics', 'off',...
                'GradObj', 'on',...
                'Display', 'off',...
                'TolFun',  1e-15,...
                'MaxFunEvals', 900);
par0=log([a;b]);            
[params, ~]=fminunc(@(params)F(params,gamma, delta,N),par0,opts);
lg_a=params(1); lg_b=params(2);                                                     
a=exp(params(1)); b=exp(lg_b); 
log_one_sigma=lg_a-lg_b;
log_one_sigma2=log_one_sigma+log(a+1)-lg_b;

function [logF, derivative]=F(params,gamma, delta, N)
a=exp(params(1)); b=exp(params(2)); %need and b to be positive 
logF=(a-N)*(log(b)-psi(a))+(b-gamma)*a/b-delta*a*(a+1)/b^2 ...
                                                      -a*log(b)+gammaln(a);
logF=-logF;                                                  
partial_a=(N-a)*psi(1,a)-gamma/b-delta*(2*a+1)/b^2+1;                                                  
ln_partial_b=-N+gamma*a/b+2*delta*a*(a+1)/b^2;   
derivative=-[partial_a*a; ln_partial_b]; 
a_=1;

% function [a,b,log_one_sigma,log_one_sigma2]=q_sigma...
%                                      (y,x,mu,l_one_w,alpha,a,b,gamma,delta)
% one_w=exp(l_one_w);                                   
% [U,D,~]=svd((bsxfun(@times,one_w,x)'*x)); d=(alpha*(1-alpha)/2)*diag(D);
% sigmaInv_d=d*a*(a+1)/b^2+1/1000;
% opts = optimset('MaxIter',1000,...
%                 'Diagnostics', 'off',...
%                 'GradObj', 'on',...
%                 'Display', 'off',...
%                 'TolFun',  1e-15,...
%                 'MaxFunEvals', 900);
% par0=log([a;b]);            
% [params, ~]=fminunc(@(params)F(params,y,x,alpha,mu,sigmaInv_d,U,d,...
%                                                   gamma, delta),par0,opts);
% % lg_a=params(1);
% lg_b=params(2);                                                     
% a=exp(params(1)); lg_a=log(a); b=exp(lg_b); 
% log_one_sigma=lg_a-lg_b;
% log_one_sigma2=log_one_sigma+log(a+1)-lg_b;
% 
% function [logF, derivative]=F(params,y,x,alpha,mu,sigmaInv_d,U,d,gamma,delta)
% a=exp(params(1)); b=exp(params(2)); %need and b to be positive
% N=length(y); 
% % sigmaInv_d=d*a*(a+1)/b^2+1/1000;
% traceTerm=x*diagMult(U,1./sqrt(sigmaInv_d),'r');
% C=sqrt((1-2*alpha)^2/4+alpha*(1-alpha))*sum(sqrt((y-x*mu).^2+sum(traceTerm.^2,2)));
% 
% %log F function (to maximise)
% logF=-0.5*sum(log((a*(a+1)/b^2)*d+1/1000))+gammaln(a)-a*psi(a)+a-sqrt(a*(a+1))/b...
%     *C-N*(log(b)-psi(a))-gamma*a/b-delta*a*(a+1)/b;
% logF=-logF; %minimise this
% 
% %derivatives wrt parameters
% partial_a=-0.5*sum((2*a+1)*d./(a*(a+1)*d+b^2/1000))-0.5*(2*a+1)*C/(b*sqrt(a*(a+1)))...
%     +(N-a)*psi(1,a)-gamma/b-delta*(2*a+1)/b+1;
% 
% partial_b=sum(a*(a+1)*d./(a*(a+1)*b*d+b^3/1000))+sqrt(a*(a+1))*C/b^2-N/b...
%     +gamma*a/b^2+delta*a*(a+1)/b^2;
% 
% %negative derivative and, multiply by a and b to account for log transform
% derivative=-[partial_a*a; partial_b*b]; 

function [alpha_,beta,l_one_w]=q_w(y,x,mu,L,l_one_sigma2,alpha)
alpha_=(1-2*alpha)^2/(2*alpha*(1-alpha))+2;
ln_beta=l_one_sigma2+log(alpha)+log((1-alpha)/2)+...
                                          log((y-x*mu).^2+sum((x*L).^2,2));
beta=exp(ln_beta);
l_one_w=0.5*(log(alpha_)-ln_beta);


function [mu, L,gamma,delta]=q_beta(y,x,l_one_sigma,l_one_sigma2,l_one_w,alpha)
% c1=log(alpha)+log((1-alpha)/2)+l_one_sigma2; 
one_w=exp(l_one_w);
% SigmaInv=exp(c1)*(bsxfun(@times,one_w,x)'*x);

%add a large variance to the diagonal- thanks to prior
% step=length(SigmaInv);
% SigmaInv(1:(step+1):end)=SigmaInv(1:(step+1):end)+1/1000;
%take the cholesky decompoistion
% L=chol(SigmaInv);
%calculate approximate posterior mean
[U,D,~]=svd((bsxfun(@times,one_w,x)'*x)); d=(alpha*(1-alpha)/2)*diag(D);
sigmaInv_d=d*exp(l_one_sigma2)+1/1000; L=diagMult(U,1./sqrt(sigmaInv_d),'r');
c1=exp(l_one_sigma2+l_one_w+log(alpha)+log((1-alpha)/2)); 
if (alpha>0.5)
    c2=-exp(log((2*alpha-1)/2)+l_one_sigma);
else
    c2=exp(log((1-2*alpha)/2)+l_one_sigma);
end
mu=c1.*y-c2;
mu=sum(bsxfun(@times,mu,x))'; %sum over all N
% mu=L\(L'\mu); %Sigma*rest_of_terms= SigmaInv^(-1)*mu
mu=U*diagMult(U',1./sigmaInv_d,'l')*mu;
gamma=-(1-2*alpha)/2*sum((y-x*mu));
delta=alpha*(1-alpha)/4*sum(one_w.*((y-x*mu).^2+sum((x*L).^2,2)));

% function [gamma_, delta,log_one_sigma,log_one_sigma2]=q_sigma...
%                                                    (y,x,mu,L,one_w,alpha,N)
% gamma_=-((1-2*alpha)/2)*sum(y-x*mu);
% delta=(alpha*(1-alpha))/4*sum(one_w.*((y-x*mu).^2+sum((x/L).^2,2)));
% log_one_sigma=log(N)-0.5*log(2*delta)+U(N+0.5,gamma_/sqrt(2*delta))-...
%                                             U(N-0.5,gamma_/sqrt(2*delta));
% if isnan(log_one_sigma), log_one_sigma=-1e6; end %give a really small value
% log_one_sigma2=log(N*(N+1))-log(2*delta)+U(N+1.5,gamma_/sqrt(2*delta))-...
%                                             U(N-0.5,gamma_/sqrt(2*delta));
% if isnan(log_one_sigma2)
%     log_one_sigma2=-1e6; end %give a really small value

% %conversion from whittakers to U
% function logU=U(a,x)
% if abs(x/sqrt(a))<1e-6
%     logU=0.5*log(pi)-sqrt(a)*x-(0.5*a+0.25)*log(2)-gammaln(0.75+0.5*a);
% elseif x/a>1e2
%     logU=-0.25*x^2-(a-0.5)*log(x);    
% else
%     logU=log(Dv(-a-0.5,x));
% end
