function [puc,lruc,pcc,lrcc,pind,lrind]=uccctests(hits,a,lag)
%[puc,tuc,pcc1,tcc1,pind1,tind1]
n=length(hits);x=sum(hits);
%LR_uc

ahat=x/n;
if x==0
   lruc=-2*n*log(1-a);
else   
   lruc=2*(x*(log(ahat)-log(a))+(n-x)*(log(1-ahat)-log(1-a)));
end
puc=1-chi2cdf(lruc,1);

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

p01=t01/(t00+t01);

if t11==0
   p11=0;
else   
   p11=t11/(t10+t11);
end

p1=(t01+t11)/n;

ll1=t00*log(1-p01);
if p01>0
    ll1=ll1+t01*log(p01)+t10*log(1-p11);
end    
if p11>0
    ll1=ll1+t11*log(p11);
end
if p1>0
   ll0=(t10+t00)*log(1-p1)+(t01+t11)*log(p1);
else
   ll0=0;
end   

lrind=2*(ll1-ll0);
pind=1-chi2cdf(lrind,1);
lrcc=lruc+lrind;
pcc=1-chi2cdf(lrcc,2);

