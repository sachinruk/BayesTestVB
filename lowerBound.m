function lowerL=lowerBound(y,mu,LSigInv,alpha,beta,a,b,quantile)
%computes the lower bound of the log likelihood
% [U, D,~]=svd(L); d=diag(D);
D=length(mu); N=length(y);
one_w=sqrt(alpha./beta); 
w=sqrt(beta./alpha)+1./alpha;
one_sigma=a/b;
muT_L=LSigInv'*mu;
muTSigmaInvmu=muT_L'*muT_L;
q1=quantile*(1-quantile);
q2=(1-2*quantile);
y2Term=q1*a*(a+1)/(2*b^2)*sum(one_w.*y.^2);

EwTerm=(q2^2/(2*q1)+2)*sum(w);

Eone_sigmaTerm=q2*one_sigma*sum(y);

ElnsigmaTerm=2*(N+1)*(log(b)-psi(a));

logData=-0.5*(N*log(2/q1)+ElnsigmaTerm+D*log(1000)...
    -muTSigmaInvmu+y2Term-Eone_sigmaTerm+EwTerm);

logQ=sum(log(diag(LSigInv)))+0.5*N*log(alpha)-N/2-gammaln(a)-log(b)+a*psi(a)...
    +psi(a)-a;

lowerL=logData-logQ;