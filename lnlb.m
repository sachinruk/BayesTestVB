function LnLik=lnlb(y,x,mu,Linv,alpha,beta,a,b,quantile)

N=length(y); D=length(mu);
q1=quantile*(1-quantile);
q2=(1-2*quantile);
one_w=sqrt(alpha./beta); w=sqrt(beta./alpha)+1./alpha;
one_sigma=a/b;
one_sigma2=a*(a+1)/b^2;
one_w_one_sigma2=one_w.*one_sigma2;
%calculate approximate posterior mean
SigmaInv=0.5*q1*bsxfun(@times,one_w_one_sigma2,x)'*x; 
%add a large variance to the diagonal- thanks to prior
SigmaInv(1:(D+1):end)=SigmaInv(1:(D+1):end)+1/1000;
L=chol(SigmaInv)';
LinvTmu=L'*mu; %L is already Sigma inv

% traceTerm=SigmaInv.*Sigma;
% traceTerm=sum(traceTerm(:));
traceTerm=Linv\L;
traceTerm=trace(traceTerm*traceTerm');
EBTSigmaInvB=LinvTmu'*LinvTmu+traceTerm;

% c1=exp(l_one_w+l_one_sigma2+log(alpha)+log((1-alpha)/2)); 
c1=one_w_one_sigma2*q1*0.5;
c2=0.5*q2*one_sigma;
muTSigmaInv=c1.*y-c2;
muTSigmaInv=sum(bsxfun(@times,muTSigmaInv,x)); %sum over all N
EmuTSigmaInvB=muTSigmaInv*mu;

y2Term=q1*a*(a+1)/(2*b^2)*sum(one_w.*y.^2);

EwTerm=(q2^2/(2*q1)+2)*sum(w);

Eone_sigmaTerm=q2*one_sigma*sum(y);

ElnsigmaTerm=2*(N+1)*(log(b)-psi(a));

logData=-0.5*(N*log(2/q1)+ElnsigmaTerm+D*log(1000)+EBTSigmaInvB+...
    -2*EmuTSigmaInvB+y2Term-Eone_sigmaTerm+EwTerm);

logQ=sum(log(diag(Linv)))-0.5*(D)+0.5*N*log(alpha)...
    -N/2-gammaln(a)-log(b)+a*psi(a)+psi(a)-a;

LnLik=logData-logQ;
