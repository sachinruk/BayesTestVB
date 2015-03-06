function mu=vb(y,x,alpha)
%uses the whittaker U function- therefore unstable

%first beta estimate
one_w=1; N=length(y);
mu= x\y; 
SigmaInv0=(alpha*(1-alpha)/2)*(x'*x); step=length(SigmaInv0);
SigmaInv0(1:(step+1):end)=SigmaInv0(1:(step+1):end)+1/1000;
L=chol(SigmaInv0);
for i=1:100
    [~,~,l_one_sigma,l_one_sigma2]=q_sigma(y,x,mu,L,one_w,alpha,N);
    [~,~,one_w]=q_w(y,x,mu,L,l_one_sigma2,alpha);
    [mu, L]=q_beta(y,x,l_one_sigma,l_one_sigma2,one_w,alpha);
end

function [mu, L]=q_beta(y,x,l_one_sigma,l_one_sigma2,l_one_w,alpha)
c1=log(alpha)+log((1-alpha)/2)+l_one_sigma2;
SigmaInv=exp(c1)*(bsxfun(@times,l_one_w,x)'*x);
%add a large variance to the diagonal- thanks to prior
step=length(SigmaInv);
SigmaInv(1:(step+1):end)=SigmaInv(1:(step+1):end)+1/1000;
%take the cholesky decompoistion
L=chol(SigmaInv);
%calculate approximate posterior mean
c1=exp(l_one_sigma2+l_one_w+log(alpha)+log((1-alpha)/2)); 
c2=exp(log((1-2*alpha)/2)+l_one_sigma);
mu=c1.*y-c2;
mu=sum(bsxfun(@times,mu,x))'; %sum over all N
mu=L\(L'\mu); %Sigma*rest_of_terms= SigmaInv^(-1)*mu

function [gamma_, delta,log_one_sigma,log_one_sigma2]=q_sigma...
                                                   (y,x,mu,L,one_w,alpha,N)
gamma_=-((1-2*alpha)/2)*sum(y-x*mu);
delta=(alpha*(1-alpha))/4*sum(one_w.*((y-x*mu).^2+sum((x/L).^2,2)));
log_one_sigma=log(N)-0.5*log(2*delta)+U(N+0.5,gamma_/sqrt(2*delta))-...
                                            U(N-0.5,gamma_/sqrt(2*delta));
if isnan(log_one_sigma), log_one_sigma=-1e6; end %give a really small value
log_one_sigma2=log(N*(N+1))-log(2*delta)+U(N+1.5,gamma_/sqrt(2*delta))-...
                                            U(N-0.5,gamma_/sqrt(2*delta));
if isnan(log_one_sigma2)
    log_one_sigma2=-1e6; end %give a really small value

function [alpha_,ln_beta,l_one_w]=q_w(y,x,mu,L,l_one_sigma2,alpha)
alpha_=(1-2*alpha)^2/(2*alpha*(1-alpha))+2;
ln_beta=l_one_sigma2+log(alpha)+log((1-alpha)/2)+...
                                          log((y-x*mu).^2+sum((x/L).^2,2));
l_one_w=0.5*(log(alpha_)-ln_beta);

%conversion from whittakers to U
function logU=U(a,x)
if abs(x/sqrt(a))<1e-6
    logU=0.5*log(pi)-sqrt(a)*x-(0.5*a+0.25)*log(2)-gammaln(0.75+0.5*a);
elseif x/a>1e2
    logU=-0.25*x^2-(a-0.5)*log(x);    
else
    logU=log(Dv(-a-0.5,x));
end