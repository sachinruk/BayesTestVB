function [zeta,pVQR] = VQRtest(returns,VaR,p,gd)
% Code to evaluate Value-at-Risk Models via Quantile Regression, from
% Gagliagone, Lima, Linton and Smith (2011)
%   'Returns' are the log returns of the data, 'VaR' are the VaR estimates
%   for a given quantile, and 'p' is the quantile of interest.

%[beta tstats VCboot itrat PseudoR2 betaboot]=quantilereg(returns,VaR,p);
[beta1 tstats itrat PseudoR2]=quantilereg(returns,VaR,p);
a0=beta1(1);a1=beta1(2);
theta=[a0 a1-1]';
T=length(returns);
muv=mean(VaR);
muv2=VaR'*VaR/T;

J=[1 muv;muv muv2];

if gd==1
   rm = RiskMetrics(returns,0.94,'Univariate');
   x=returns./rm;
   q5=a0+a1*VaR;q2=q5./rm;
   f=ksdensity(x,q2);
   fx=f./rm;
elseif gd==2
   %fit a GARCH(1,1)-N model
   spec1 = garchset('P',1,'Q',1,'Display','Off');
   [cGC,errors,LLF,innovationsC,rm] = garchfit(spec1,returns); 
   x=returns./rm;
   q5=a0+a1*VaR;q2=q5./rm;
   f=ksdensity(x,q2);
   fx=f./rm;
else
    q5=a0+a1*VaR;
   % Koenker and MAchado (1999)
   hn=T^(-1/3)*(norminv(1-p/2)^(2))^(1/3)*(1.5*normpdf(norminv(p))^2/(2*norminv(p)^2+1))^(1/3);
   if hn>p
      hn=p/2;
   end   
   % beta(p+hn), beta(p-hn)
   [betap tstats1 itrat1 PseudoR21]=quantilereg(returns,VaR,p+hn);
   [betam tstats1 itrat1 PseudoR21]=quantilereg(returns,VaR,p-hn);
   q5p=betap(1)+betap(2)*VaR;q5m=betam(1)+betam(2)*VaR;
   di=abs(q5p-q5m);di(di==0)=min(di(di>0));
   sf=di/(2*hn);fx=1./sf;  %conditional density (on VaR), assuming iid errors
%    del=std(returns)/mean(sf);
%    rm=sf/del;rm1=rm./mean(rm);rm2=rm1./std(rm1)*(std(returns));%+std(returns);
%    st1=1;
%    if abs(mean(rm2)-std(returns))>0.05
%        st1=0;
%    end    
%    while st1==0
%        if abs(mean(rm2)-std(returns))<0.05
%           st1=1;
%        else
%           if mean(rm2)>std(returns)
%              rm2=rm2./1.5;
%           else
%              rm2=rm2*1.2;
%           end 
%        end   
%    end 
%    %rm2=rm;
%    x=returns./rm2;q2=q5./rm2;
%    f=ksdensity(x,q2);fx=f./rm2;
   
end


mud=mean(fx);
mudv=(VaR'*fx)/T;
% n=length(returns);
% dv2=zeros(n,1);
% for i=1:n
%     dv2(i)=VaR(i)^2*fx(i);
% end
mudv2=mean(VaR.^2.*fx);%mean(dv2);

H=[mud mudv;mudv mudv2];

iH=H\eye(2,2);
iHJiH=iH*J*iH;

zeta=T*(theta'*(iHJiH\eye(2,2))*theta)/(p*(1-p));

%zeta=T*(theta'*inv(p*(1-p)*inv(H)*J*inv(H))*theta);
pVQR=1-chi2cdf(zeta,2);

end

