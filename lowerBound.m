function lowerL=lowerBound(y,mu,LSigInv,alpha_,beta,a,b,quantile)
%computes the lower bound of the log likelihood
% [U, D,~]=svd(L); d=diag(D);
D=length(mu); N=length(y);
one_w=sqrt(alpha_./beta); w=sqrt(beta./alpha_)+1./alpha_;
one_sigma=a/b;
muT_L=LSigInv'*mu;
muTSigmaInvmu=muT_L'*muT_L;
logData=-0.5*(N*log(2/(quantile*(1-quantile)))+2*(N+1)*(log(b)-psi(a))+D*log(1000)...
    -muTSigmaInvmu+quantile*(1-quantile)*a*(a+1)/(2*b^2)*sum(one_w.*y.^2)...
    +((1-2*quantile)^2/(2*quantile*(1-quantile))+2)*sum(w)-(1-2*quantile)*one_sigma*sum(y));

logQ=sum(log(diag(LSigInv)))+0.5*sum(log(alpha_))-N/2-gammaln(a)-log(b)+a*psi(a)...
    +psi(a)-a;

lowerL=logData-logQ;