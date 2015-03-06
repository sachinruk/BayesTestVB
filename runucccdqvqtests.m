function [puc,tuc,pcc1,tcc1,pind1,tind1,pcc4,tcc4,pind4,tind4,pdq1,dq1,pdq4,dq4,tsvq,pvqr]=runucccdqvqtests(y,f,p)

gd=3;
hits=(y<f);
[puc,tuc,pcc1,tcc1,pind1,tind1]=uccctests(hits,p,1);
[puc,tuc,pcc4,tcc4,pind4,tind4]=uccctests(hits,p,4);
[pdq1,dq1]=dqtest(y',f',p,1);[pdq4,dq4]=dqtest(y',f',p,4);
[tsvq,pvqr] = VQRtest(y',f',p,gd);