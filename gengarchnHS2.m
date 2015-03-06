function [y, sig, VaR5, VaR1, VaR5HS, VaR1HS]=gengarchnHS2(n,a0,a1,b1)

l=500;
sig=zeros(1,n+l); 
% VaR1=zeros(1,n+l); VaR5=zeros(1,n+l); 
VaR1HS=zeros(1,n+l); VaR5HS=zeros(1,n+l);
y=randn(1,n+l);
%y(1)=randn*sig;
for t=2:n+l
    sig(t)=sqrt(a0+a1*y(t-1)^2+b1*sig(t-1)^2);
    y(t)=y(t)*sig(t);
    if t>500
        p_tiles=prctile(y(t-250:t-1),[1 5]);
        VaR1HS(t)=p_tiles(1); VaR5HS(t)=p_tiles(2);
%         VaR5HS(t)=prctile(y(t-250:t-1),5);
    end
end

%the following does not need to be in a loop
VaR1=sig*norminv(0.01);
VaR5=sig*norminv(0.05);

%burnouts
sig=sig(501:end);
y=y(501:end);
VaR1=VaR1(501:end);VaR1HS=VaR1HS(501:end);
VaR5=VaR5(501:end);VaR5HS=VaR5HS(501:end);