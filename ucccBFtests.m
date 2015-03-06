function [BFuc,BFcc,BFind]=ucccBFtests(hits,a,lag)

n=length(hits);x=sum(hits);
%BF_uc

lbf=gammaln(n-x+1)+gammaln(x+1)-gammaln(n+2)-x*log(a)-(n-x)*log(1-a);
BFuc = exp(lbf);

%% Ind test
r5=hits(lag+1:n);r51=0*r5;
for k=1:n-lag
    r51(k)=max(hits(k:k+lag-1));
end
i11=r5.*r51;%(r5<V5(2:500)').*(r5l<V5(1:499)');
i01=r5.*(1-r51);%(r5<V5(2:500)').*(r5l>V5(1:499)');
i10=(1-r5).*r51;%(r5>V5(2:500)').*(r5l<V5(1:499)');
i00=(1-r5).*(1-r51);%(r5>V5(2:500)').*(r5l>V5(1:499)');

t00=sum(i00);t01=sum(i01);t10=sum(i10);t11=sum(i11);

n1=t00+t10+t11+t01;
lbf1=gammaln(t01+1)+gammaln(t00+1)-gammaln(t01+t00+2)+gammaln(t11+1)+gammaln(t10+1)-gammaln(t11+t10+2)-(gammaln(t00+t10+1)+gammaln(t11+t01+1)-gammaln(n1+2));
BFind = exp(lbf1);

x1=t01+t11;x2=t10+t00;
lbf2=gammaln(t00+1)+gammaln(t01+1)-gammaln(t01+t00+2)+gammaln(t10+1)+gammaln(t11+1)-gammaln(t10+t11+2)-x1*log(a)-x2*log(1-a);
BFcc = exp(lbf2);

