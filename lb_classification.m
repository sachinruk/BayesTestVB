function l=lb_classification(X,t,mu,SigmaInv,C)

N=length(t);
d=length(mu);
[U, D,~]=svd(SigmaInv); d_=diag(D); eigSigma=1./d_;
TrSigma=sum(eigSigma);
TrXSigmaX=diagMult(X*U,sqrt(eigSigma),'r'); TrXSigmaX=sum(TrXSigmaX(:).^2);
lnP=-N/2*log(2*pi)-d/2*log(1000)-1/(2*1000)*(mu'*mu+TrSigma)-0.5*TrXSigmaX;
lnBeta=-0.5*sum(log(eigSigma))-d/2;
lny=-sum(log(C));

l=lnP-lnBeta-lny;