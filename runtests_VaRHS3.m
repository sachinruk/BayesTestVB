% generate data

a0=0.1;a1=0.1;b1=0.85;

n=2500
for j=1:5000
disp(j);
[y, sig, VaR5, VaR1, VaR5HS, VaR1HS]=gengarchnHS2(n,a0,a1,b1);
  hits5s2500(j)=sum(y<VaR5);hits5p2500(j)=sum(y<VaR5HS);
  hits1s2500(j)=sum(y<VaR1);hits1p2500(j)=sum(y<VaR1HS);
  [pucn1s2500(j),lrucn1s2500(j),pccn1s2500(j),lrccn1s2500(j),pindn1s2500(j),lrindn1s2500(j),pcc4n1s2500(j),lrcc4n1s2500(j),pind4n1s2500(j),lrind4n1s2500(j),pdq1n1s2500(j),dq1n1s2500(j),pdq4n1s2500(j),dq4n1s2500(j),tsvqn1s2500(j),pvqrn1s2500(j)]=runucccdqvqtests(y,VaR1,0.01);
  [pucn5s2500(j),lrucn5s2500(j),pccn5s2500(j),lrccn5s2500(j),pindn5s2500(j),lrindn5s2500(j),pcc4n5s2500(j),lrcc4n5s2500(j),pind4n5s2500(j),lrind4n5s2500(j),pdq1n5s2500(j),dq1n5s2500(j),pdq4n5s2500(j),dq4n5s2500(j),tsvqn5s2500(j),pvqrn5s2500(j)]=runucccdqvqtests(y,VaR5,0.05);
  [pucn1p2500(j),lrucn1p2500(j),pccn1p2500(j),lrccn1p2500(j),pindn1p2500(j),lrindn1p2500(j),pcc4n1p2500(j),lrcc4n1p2500(j),pind4n1p2500(j),lrind4n1p2500(j),pdq1n1p2500(j),dq1n1p2500(j),pdq4n1p2500(j),dq4n1p2500(j),tsvqn1p2500(j),pvqrn1p2500(j)]=runucccdqvqtests(y,VaR1HS,0.01);
  [pucn5p2500(j),lrucn5p2500(j),pccn5p2500(j),lrccn5p2500(j),pindn5p2500(j),lrindn5p2500(j),pcc4n5p2500(j),lrcc4n5p2500(j),pind4n5p2500(j),lrind4n5p2500(j),pdq1n5p2500(j),dq1n5p2500(j),pdq4n5p2500(j),dq4n5p2500(j),tsvqn5p2500(j),pvqrn5p2500(j)]=runucccdqvqtests(y,VaR5HS,0.05);
  [BFucn1s2500(j),BFcc1n1s2500(j),BFind1n1s2500(j),BFcc4n1s2500(j),BFind4n1s2500(j),BFvqn1s2500(j)]=runBFtests(y,VaR1,0.01);
  [BFucn5s2500(j),BFcc1n5s2500(j),BFind1n5s2500(j),BFcc4n5s2500(j),BFind4n5s2500(j),BFvqn5s2500(j)]=runBFtests(y,VaR5,0.05);
 [BFucn1p2500(j),BFcc1n1p2500(j),BFind1n1p2500(j),BFcc4n1p2500(j),BFind4n1p2500(j),BFvqn1p2500(j)]=runBFtests(y,VaR1HS,0.01);
 [BFucn5p2500(j),BFcc1n5p2500(j),BFind1n5p2500(j),BFcc4n5p2500(j),BFind4n5p2500(j),BFvqn5p2500(j)]=runBFtests(y,VaR5HS,0.05);
  BFdq1n1s2500(j) = 1/bfdqmlik3(y,VaR1,0.01);BFdq1n5s2500(j)=1/bfdqmlik3(y,VaR5,0.05);
  BFdq1n1p2500(j) = 1/bfdqmlik3(y,VaR1HS,0.01);BFdq1n5p2500(j)=1/bfdqmlik3(y,VaR5HS,0.05);
  BFdq14n1s2500(j) = 1/bfdqmlik4(y,VaR1,0.01);BFdq14n5s2500(j)=1/bfdqmlik4(y,VaR5,0.05);
  BFdq14n1p2500(j) = 1/bfdqmlik4(y,VaR1HS,0.01);BFdq14n5p2500(j)=1/bfdqmlik4(y,VaR5HS,0.05);
   BFdq2n1s2500(j) = 1/bfdq2mlik2(y,VaR1,0.01);BFdq2n5s2500(j)=1/bfdq2mlik2(y,VaR5,0.05);
   BFdq2n1p2500(j) = 1/bfdq2mlik2(y,VaR1HS,0.01);BFdq2n5p2500(j)=1/bfdq2mlik2(y,VaR5HS,0.05);
  
  if rem(j,100)==0
     disp(j);
     [pucn1s2500(j),lrucn1s2500(j),pccn1s2500(j),lrccn1s2500(j),pindn1s2500(j),lrindn1s2500(j),pcc4n1s2500(j),lrcc4n1s2500(j),pind4n1s2500(j),lrind4n1s2500(j),pdq1n1s2500(j),dq1n1s2500(j),pdq4n1s2500(j),dq4n1s2500(j),tsvqn1s2500(j),pvqrn1s2500(j);
      pucn1p2500(j),lrucn1p2500(j),pccn1p2500(j),lrccn1p2500(j),pindn1p2500(j),lrindn1p2500(j),pcc4n1p2500(j),lrcc4n1p2500(j),pind4n1p2500(j),lrind4n1p2500(j),pdq1n1p2500(j),dq1n1p2500(j),pdq4n1p2500(j),dq4n1p2500(j),tsvqn1p2500(j),pvqrn1p2500(j)]
     %[2*log(BFdq2n1s2500(j)) 2*log(BFdq2n5s2500(j)) 2*log(BFdq14n1s2500(j)) 2*log(BFdq14n5s2500(j));
     %  2*log(BFdq2n1p2500(j)) 2*log(BFdq2n5p2500(j)) 2*log(BFdq14n1p2500(j)) 2*log(BFdq14n5p2500(j))]  
  end   
  
end  

disp('size');
[length(pucn5s2500(pucn5s2500<0.05)) length(pindn5s2500(pindn5s2500<0.05)) length(pccn5s2500(pccn5s2500<0.05)) length(pind4n5s2500(pind4n5s2500<0.05)) length(pcc4n5s2500(pcc4n5s2500<0.05)) length(pdq1n5s2500(pdq1n5s2500<0.05)) length(pdq4n5s2500(pdq4n5s2500<0.05)) length(pvqrn5s2500(pvqrn5s2500<0.05)) length(BFucn5s2500(BFucn5s2500>1)) length(BFind1n5s2500(BFind1n5s2500>1)) length(BFcc1n5s2500(BFcc1n5s2500>1)) length(BFind4n5s2500(BFind4n5s2500>1)) length(BFcc4n5s2500(BFcc4n5s2500>1)) length(BFvqn5s2500(BFvqn5s2500>1)) length(BFdq1n5s2500(BFdq1n5s2500>1)) length(BFdq14n5s2500(BFdq14n5s2500>1)) length(BFdq2n5s2500(BFdq2n5s2500>1))]./length(pvqrn5s2500)
[length(pucn1s2500(pucn1s2500<0.05)) length(pindn1s2500(pindn1s2500<0.05)) length(pccn1s2500(pccn1s2500<0.05)) length(pind4n1s2500(pind4n1s2500<0.05)) length(pcc4n1s2500(pcc4n1s2500<0.05)) length(pdq1n1s2500(pdq1n1s2500<0.05)) length(pdq4n1s2500(pdq4n1s2500<0.05)) length(pvqrn1s2500(pvqrn1s2500<0.05)) length(BFucn1s2500(BFucn1s2500>1)) length(BFind1n1s2500(BFind1n1s2500>1)) length(BFcc1n1s2500(BFcc1n1s2500>1))  length(BFind4n1s2500(BFind4n1s2500>1)) length(BFcc4n1s2500(BFcc4n1s2500>1)) length(BFvqn1s2500(BFvqn1s2500>1)) length(BFdq1n1s2500(BFdq1n1s2500>1)) length(BFdq14n1s2500(BFdq14n1s2500>1)) length(BFdq2n1s2500(BFdq2n1s2500>1))]./length(pvqrn1s2500)
%[length(pucn5s2500(pucn5s2500<0.05)) length(pindn5s2500(pindn5s2500<0.05)) length(pccn5s2500(pccn5s2500<0.05)) length(pind4n5s2500(pind4n5s2500<0.05)) length(pcc4n5s2500(pcc4n5s2500<0.05)) length(pdq1n5s2500(pdq1n5s2500<0.05)) length(pdq4n5s2500(pdq4n5s2500<0.05)) length(pvqrn5s2500(pvqrn5s2500<0.05)) length(BFucn5s2500(BFucn5s2500>1)) length(BFind1n5s2500(BFind1n5s2500>1)) length(BFcc1n5s2500(BFcc1n5s2500>1)) length(BFind4n5s2500(BFind4n5s2500>1)) length(BFcc4n5s2500(BFcc4n5s2500>1)) length(BFvqn5s2500(BFvqn5s2500>1)) length(BFdq1n5s2500(BFdq1n5s2500>1)) length(BFdq14n5s2500(BFdq14n5s2500>1))]./length(pvqrn5s2500)
%[length(pucn1s2500(pucn1s2500<0.05)) length(pindn1s2500(pindn1s2500<0.05)) length(pccn1s2500(pccn1s2500<0.05)) length(pind4n1s2500(pind4n1s2500<0.05)) length(pcc4n1s2500(pcc4n1s2500<0.05)) length(pdq1n1s2500(pdq1n1s2500<0.05)) length(pdq4n1s2500(pdq4n1s2500<0.05)) length(pvqrn1s2500(pvqrn1s2500<0.05)) length(BFucn1s2500(BFucn1s2500>1)) length(BFind1n1s2500(BFind1n1s2500>1)) length(BFcc1n1s2500(BFcc1n1s2500>1))  length(BFind4n1s2500(BFind4n1s2500>1)) length(BFcc4n1s2500(BFcc4n1s2500>1)) length(BFvqn1s2500(BFvqn1s2500>1)) length(BFdq1n1s2500(BFdq1n1s2500>1)) length(BFdq14n1s2500(BFdq14n1s2500>1))]./length(pvqrn1s2500)


disp('power');
[length(pucn5p2500(pucn5p2500<0.05)) length(pindn5p2500(pindn5p2500<0.05)) length(pccn5p2500(pccn5p2500<0.05)) length(pind4n5p2500(pind4n5p2500<0.05)) length(pcc4n5p2500(pcc4n5p2500<0.05)) length(pdq1n5p2500(pdq1n5p2500<0.05)) length(pdq4n5p2500(pdq4n5p2500<0.05)) length(pvqrn5p2500(pvqrn5p2500<0.05)) length(BFucn5p2500(BFucn5p2500>1)) length(BFind1n5p2500(BFind1n5p2500>1)) length(BFcc1n5p2500(BFcc1n5p2500>1)) length(BFind4n5p2500(BFind4n5p2500>1)) length(BFcc4n5p2500(BFcc4n5p2500>1)) length(BFvqn5p2500(BFvqn5p2500>1)) length(BFdq1n5p2500(BFdq1n5p2500>1)) length(BFdq14n5p2500(BFdq14n5p2500>1)) length(BFdq2n5p2500(BFdq2n5p2500>1))]./length(pvqrn5p2500)
[length(pucn1p2500(pucn1p2500<0.05)) length(pindn1p2500(pindn1p2500<0.05)) length(pccn1p2500(pccn1p2500<0.05)) length(pind4n1p2500(pind4n1p2500<0.05)) length(pcc4n1p2500(pcc4n1p2500<0.05)) length(pdq1n1p2500(pdq1n1p2500<0.05)) length(pdq4n1p2500(pdq4n1p2500<0.05)) length(pvqrn1p2500(pvqrn1p2500<0.05)) length(BFucn1p2500(BFucn1p2500>1)) length(BFind1n1p2500(BFind1n1p2500>1)) length(BFcc1n1p2500(BFcc1n1p2500>1)) length(BFind4n1p2500(BFind4n1p2500>1)) length(BFcc4n1p2500(BFcc4n1p2500>1)) length(BFvqn1p2500(BFvqn1p2500>1)) length(BFdq1n1p2500(BFdq1n1p2500>1)) length(BFdq14n1p2500(BFdq14n5p2500>1)) length(BFdq2n1p2500(BFdq2n1p2500>1))]./length(pvqrn1p2500)
%[length(pucn5p2500(pucn5p2500<0.05)) length(pindn5p2500(pindn5p2500<0.05)) length(pccn5p2500(pccn5p2500<0.05)) length(pind4n5p2500(pind4n5p2500<0.05)) length(pcc4n5p2500(pcc4n5p2500<0.05)) length(pdq1n5p2500(pdq1n5p2500<0.05)) length(pdq4n5p2500(pdq4n5p2500<0.05)) length(pvqrn5p2500(pvqrn5p2500<0.05)) length(BFucn5p2500(BFucn5p2500>1)) length(BFind1n5p2500(BFind1n5p2500>1)) length(BFcc1n5p2500(BFcc1n5p2500>1)) length(BFind4n5p2500(BFind4n5p2500>1)) length(BFcc4n5p2500(BFcc4n5p2500>1)) length(BFvqn5p2500(BFvqn5p2500>1)) length(BFdq1n5p2500(BFdq1n5p2500>1)) length(BFdq14n5p2500(BFdq14n5p2500>1))]./length(pvqrn5p2500)
%[length(pucn1p2500(pucn1p2500<0.05)) length(pindn1p2500(pindn1p2500<0.05)) length(pccn1p2500(pccn1p2500<0.05)) length(pind4n1p2500(pind4n1p2500<0.05)) length(pcc4n1p2500(pcc4n1p2500<0.05)) length(pdq1n1p2500(pdq1n1p2500<0.05)) length(pdq4n1p2500(pdq4n1p2500<0.05)) length(pvqrn1p2500(pvqrn1p2500<0.05)) length(BFucn1p2500(BFucn1p2500>1)) length(BFind1n1p2500(BFind1n1p2500>1)) length(BFcc1n1p2500(BFcc1n1p2500>1)) length(BFind4n1p2500(BFind4n1p2500>1)) length(BFcc4n1p2500(BFcc4n1p2500>1)) length(BFvqn1p2500(BFvqn1p2500>1)) length(BFdq1n1p2500(BFdq1n1p2500>1)) length(BFdq14n1p2500(BFdq14n1p2500>1))]./length(pvqrn1p2500)

BFuc5s2500c=prctile(BFucn5s2500,95);BFind5s2500c=prctile(BFind1n5s2500,95);
BFcc5s2500c=prctile(BFcc1n5s2500,95);BFcc45s2500c=prctile(BFcc4n5s2500,95);
BFind45s2500c=prctile(BFind4n5s2500,95);BFvq5s2500c=prctile(BFvqn5s2500,95);
BFdq15s2500c=prctile(BFdq1n5s2500,95);BFdq145s2500c=prctile(BFdq14n5s2500,95);BFdq25s2500c=prctile(BFdq2n5s2500,95);
uc5s2500c=prctile(pucn5s2500,5);ind5s2500c=prctile(pindn5s2500,5);ind45s2500c=prctile(pind4n5s2500,5);
cc5s2500c=prctile(pccn5s2500,5);cc45s2500c=prctile(pcc4n5s2500,5);
dq15s2500c=prctile(dq1n5s2500,95);dq45s2500c=prctile(dq4n5s2500,95);tsvq5s2500c=prctile(tsvqn5s2500,95);

BFuc1s2500c=prctile(BFucn1s2500,95);BFind1s2500c=prctile(BFind1n1s2500,95);
BFcc1s2500c=prctile(BFcc1n1s2500,95);BFcc41s2500c=prctile(BFcc4n1s2500,95);
BFind41s2500c=prctile(BFind4n1s2500,95);BFvq1s2500c=prctile(BFvqn1s2500,95);
BFdq11s2500c=prctile(BFdq1n1s2500,95);BFdq141s2500c=prctile(BFdq14n1s2500,95);BFdq21s2500c=prctile(BFdq2n1s2500,95);
uc1s2500c=prctile(pucn1s2500,5);ind1s2500c=prctile(pindn1s2500,5);ind41s2500c=prctile(pind4n1s2500,5);
cc1s2500c=prctile(pccn1s2500,5);cc41s2500c=prctile(pcc4n1s2500,5);
dq11s2500c=prctile(dq1n1s2500,95);dq41s2500c=prctile(dq4n1s2500,95);tsvq1s2500c=prctile(tsvqn1s2500,95);

disp('corrected size');
[length(pucn5s2500(pucn5s2500<=uc5s2500c)) length(pindn5s2500(pindn5s2500<ind5s2500c)) length(pccn5s2500(pccn5s2500<cc5s2500c)) length(pind4n5s2500(pind4n5s2500<ind45s2500c)) length(pcc4n5s2500(pcc4n5s2500<cc45s2500c)) length(dq1n5s2500(dq1n5s2500>dq15s2500c)) length(dq4n5s2500(dq4n5s2500>dq45s2500c)) length(tsvqn5s2500(tsvqn5s2500>tsvq5s2500c)) length(BFucn5s2500(BFucn5s2500>BFuc5s2500c)) length(BFind1n5s2500(BFind1n5s2500>BFind5s2500c)) length(BFcc1n5s2500(BFcc1n5s2500>BFcc5s2500c)) length(BFind4n5s2500(BFind4n5s2500>BFind45s2500c)) length(BFcc4n5s2500(BFcc4n5s2500>BFcc45s2500c)) length(BFvqn5s2500(BFvqn5s2500>BFvq5s2500c)) length(BFdq1n5s2500(BFdq1n5s2500>BFdq15s2500c)) length(BFdq14n5s2500(BFdq14n5s2500>BFdq145s2500c)) length(BFdq2n5s2500(BFdq2n5s2500>BFdq25s2500c))]./length(pvqrn5s2500)
[length(pucn1s2500(pucn1s2500<=uc1s2500c)) length(pindn1s2500(pindn1s2500<ind1s2500c)) length(pccn1s2500(pccn1s2500<cc1s2500c)) length(pind4n1s2500(pind4n1s2500<ind41s2500c)) length(pcc4n1s2500(pcc4n1s2500<cc41s2500c)) length(dq1n1s2500(dq1n1s2500>dq11s2500c)) length(dq4n1s2500(dq4n1s2500>dq41s2500c)) length(tsvqn1s2500(tsvqn1s2500>tsvq1s2500c)) length(BFucn1s2500(BFucn1s2500>BFuc1s2500c)) length(BFind1n1s2500(BFind1n1s2500>BFind1s2500c)) length(BFcc1n1s2500(BFcc1n1s2500>BFcc1s2500c)) length(BFind4n1s2500(BFind4n1s2500>BFind41s2500c)) length(BFcc4n1s2500(BFcc4n1s2500>BFcc41s2500c)) length(BFvqn1s2500(BFvqn1s2500>BFvq1s2500c)) length(BFdq1n1s2500(BFdq1n1s2500>BFdq11s2500c)) length(BFdq14n1s2500(BFdq14n1s2500>BFdq141s2500c)) length(BFdq2n1s2500(BFdq2n1s2500>BFdq21s2500c))]./length(pvqrn1s2500)
%[length(pucn5s2500(pucn5s2500<=uc5s2500c)) length(pindn5s2500(pindn5s2500<ind5s2500c)) length(pccn5s2500(pccn5s2500<cc5s2500c)) length(pind4n5s2500(pind4n5s2500<ind45s2500c)) length(pcc4n5s2500(pcc4n5s2500<cc45s2500c)) length(dq1n5s2500(dq1n5s2500>dq15s2500c)) length(dq4n5s2500(dq4n5s2500>dq45s2500c)) length(tsvqn5s2500(tsvqn5s2500>tsvq5s2500c)) length(BFucn5s2500(BFucn5s2500>BFuc5s2500c)) length(BFind1n5s2500(BFind1n5s2500>BFind5s2500c)) length(BFcc1n5s2500(BFcc1n5s2500>BFcc5s2500c)) length(BFind4n5s2500(BFind4n5s2500>BFind45s2500c)) length(BFcc4n5s2500(BFcc4n5s2500>BFcc45s2500c)) length(BFvqn5s2500(BFvqn5s2500>BFvq5s2500c)) length(BFdq1n5s2500(BFdq1n5s2500>BFdq15s2500c)) length(BFdq14n5s2500(BFdq14n5s2500>BFdq145s2500c))]./length(pvqrn5s2500)
%[length(pucn1s2500(pucn1s2500<=uc1s2500c)) length(pindn1s2500(pindn1s2500<ind1s2500c)) length(pccn1s2500(pccn1s2500<cc1s2500c)) length(pind4n1s2500(pind4n1s2500<ind41s2500c)) length(pcc4n1s2500(pcc4n1s2500<cc41s2500c)) length(dq1n1s2500(dq1n1s2500>dq11s2500c)) length(dq4n1s2500(dq4n1s2500>dq41s2500c)) length(tsvqn1s2500(tsvqn1s2500>tsvq1s2500c)) length(BFucn1s2500(BFucn1s2500>BFuc1s2500c)) length(BFind1n1s2500(BFind1n1s2500>BFind1s2500c)) length(BFcc1n1s2500(BFcc1n1s2500>BFcc1s2500c)) length(BFind4n1s2500(BFind4n1s2500>BFind41s2500c)) length(BFcc4n1s2500(BFcc4n1s2500>BFcc41s2500c)) length(BFvqn1s2500(BFvqn1s2500>BFvq1s2500c)) length(BFdq1n1s2500(BFdq1n1s2500>BFdq11s2500c)) length(BFdq14n1s2500(BFdq14n1s2500>BFdq141s2500c))]./length(pvqrn1s2500)


disp('size-adjusted power');
[length(pucn5p2500(pucn5p2500<=uc5s2500c)) length(pindn5p2500(pindn5p2500<ind5s2500c)) length(pccn5p2500(pccn5p2500<cc5s2500c)) length(pind4n5p2500(pind4n5p2500<ind45s2500c)) length(pcc4n5p2500(pcc4n5p2500<cc45s2500c)) length(dq1n5p2500(dq1n5p2500>dq15s2500c)) length(dq4n5p2500(dq4n5p2500>dq45s2500c)) length(tsvqn5p2500(tsvqn5p2500>tsvq5s2500c)) length(BFucn5p2500(BFucn5p2500>BFuc5s2500c)) length(BFind1n5p2500(BFind1n5p2500>BFind5s2500c)) length(BFcc1n5p2500(BFcc1n5p2500>BFcc5s2500c)) length(BFind4n5p2500(BFind4n5p2500>BFind45s2500c)) length(BFcc4n5p2500(BFcc4n5p2500>BFcc45s2500c)) length(BFvqn5p2500(BFvqn5p2500>BFvq5s2500c)) length(BFdq1n5p2500(BFdq1n5p2500>BFdq15s2500c)) length(BFdq14n5p2500(BFdq14n5p2500>BFdq145s2500c)) length(BFdq2n5p2500(BFdq2n5p2500>BFdq25s2500c))]./length(pvqrn5p2500)
[length(pucn1p2500(pucn1p2500<=uc1s2500c)) length(pindn1p2500(pindn1p2500<ind1s2500c)) length(pccn1p2500(pccn1p2500<cc1s2500c)) length(pind4n1p2500(pind4n1p2500<ind41s2500c)) length(pcc4n1p2500(pcc4n1p2500<cc41s2500c)) length(dq1n1p2500(dq1n1p2500>dq11s2500c)) length(dq4n1p2500(dq4n1p2500>dq41s2500c)) length(tsvqn1p2500(tsvqn1p2500>tsvq1s2500c)) length(BFucn1p2500(BFucn1p2500>BFuc1s2500c)) length(BFind1n1p2500(BFind1n1p2500>BFind1s2500c)) length(BFcc1n1p2500(BFcc1n1p2500>BFcc1s2500c)) length(BFind4n1p2500(BFind4n1p2500>BFind41s2500c)) length(BFcc4n1p2500(BFcc4n1p2500>BFcc41s2500c)) length(BFvqn1p2500(BFvqn1p2500>BFvq1s2500c)) length(BFdq1n1p2500(BFdq1n1p2500>BFdq11s2500c)) length(BFdq14n1p2500(BFdq14n1p2500>BFdq141s2500c)) length(BFdq2n1p2500(BFdq2n1p2500>BFdq21s2500c))]./length(pvqrn1p2500)
%[length(pucn5p2500(pucn5p2500<=uc5s2500c)) length(pindn5p2500(pindn5p2500<ind5s2500c)) length(pccn5p2500(pccn5p2500<cc5s2500c)) length(pind4n5p2500(pind4n5p2500<ind45s2500c)) length(pcc4n5p2500(pcc4n5p2500<cc45s2500c)) length(dq1n5p2500(dq1n5p2500>dq15s2500c)) length(dq4n5p2500(dq4n5p2500>dq45s2500c)) length(tsvqn5p2500(tsvqn5p2500>tsvq5s2500c)) length(BFucn5p2500(BFucn5p2500>BFuc5s2500c)) length(BFind1n5p2500(BFind1n5p2500>BFind5s2500c)) length(BFcc1n5p2500(BFcc1n5p2500>BFcc5s2500c)) length(BFind4n5p2500(BFind4n5p2500>BFind45s2500c)) length(BFcc4n5p2500(BFcc4n5p2500>BFcc45s2500c)) length(BFvqn5p2500(BFvqn5p2500>BFvq5s2500c)) length(BFdq1n5p2500(BFdq1n5p2500>BFdq15s2500c)) length(BFdq14n5p2500(BFdq14n5p2500>BFdq145s2500c))]./length(pvqrn5p2500)
%[length(pucn1p2500(pucn1p2500<=uc1s2500c)) length(pindn1p2500(pindn1p2500<ind1s2500c)) length(pccn1p2500(pccn1p2500<cc1s2500c)) length(pind4n1p2500(pind4n1p2500<ind41s2500c)) length(pcc4n1p2500(pcc4n1p2500<cc41s2500c)) length(dq1n1p2500(dq1n1p2500>dq11s2500c)) length(dq4n1p2500(dq4n1p2500>dq41s2500c)) length(tsvqn1p2500(tsvqn1p2500>tsvq1s2500c)) length(BFucn1p2500(BFucn1p2500>BFuc1s2500c)) length(BFind1n1p2500(BFind1n1p2500>BFind1s2500c)) length(BFcc1n1p2500(BFcc1n1p2500>BFcc1s2500c)) length(BFind4n1p2500(BFind4n1p2500>BFind41s2500c)) length(BFcc4n1p2500(BFcc4n1p2500>BFcc41s2500c)) length(BFvqn1p2500(BFvqn1p2500>BFvq1s2500c)) length(BFdq1n1p2500(BFdq1n1p2500>BFdq11s2500c)) length(BFdq14n1p2500(BFdq14n1p2500>BFdq141s2500c))]./length(pvqrn1p2500)

Ms52500=[hits5s2500' pucn5s2500' lrucn5s2500' pindn5s2500' lrindn5s2500' pccn5s2500' lrccn5s2500' pind4n5s2500' lrind4n5s2500' pcc4n5s2500' lrcc4n5s2500' dq1n5s2500' pdq1n5s2500' dq4n5s2500' pdq4n5s2500' tsvqn5s2500' pvqrn5s2500' BFucn5s2500' BFind1n5s2500' BFcc1n5s2500' BFind4n5s2500' BFcc4n5s2500' BFvqn5s2500' BFdq1n5s2500' BFdq14n5s2500'];% BFdq2n5s2500']; 
%Ms52500=[pucn5s2500' lrucn5s2500' pindn5s2500' lrindn5s2500' pccn5s2500' lrccn5s2500' pind4n5s2500' lrind4n5s2500' pcc4n5s2500' lrcc4n5s2500' dq1n5s2500' pdq1n5s2500' dq4n5s2500' pdq4n5s2500' tsvqn5s2500' pvqrn5s2500' BFucn5s2500' BFind1n5s2500' BFcc1n5s2500' BFind4n5s2500' BFcc4n5s2500' BFvqn5s2500']; 
csvwrite('HS2500s_5final3.csv',Ms52500)
Ms12500=[hits1s2500' pucn1s2500' lrucn1s2500' pindn1s2500' lrindn1s2500' pccn1s2500' lrccn1s2500' pind4n1s2500' lrind4n1s2500' pcc4n1s2500' lrcc4n1s2500' dq1n1s2500' pdq1n1s2500' dq4n1s2500' pdq4n1s2500' tsvqn1s2500' pvqrn1s2500' BFucn1s2500' BFind1n1s2500' BFcc1n1s2500' BFind4n1s2500' BFcc4n1s2500' BFvqn1s2500' BFdq1n1s2500' BFdq14n1s2500'];% BFdq2n1s2500']; 
%Ms12500=[pucn1s2500' lrucn1s2500' pindn1s2500' lrindn1s2500' pccn1s2500' lrccn1s2500' pind4n1s2500' lrind4n1s2500' pcc4n1s2500' lrcc4n1s2500' dq1n1s2500' pdq1n1s2500' dq4n1s2500' pdq4n1s2500' tsvqn1s2500' pvqrn1s2500' BFucn1s2500' BFind1n1s2500' BFcc1n1s2500' BFind4n1s2500' BFcc4n1s2500' BFvqn1s2500']; 
csvwrite('HS2500s_1final3.csv',Ms12500)
Mp52500=[hits5p2500' pucn5p2500' lrucn5p2500' pindn5p2500' lrindn5p2500' pccn5p2500' lrccn5p2500' pind4n5p2500' lrind4n5p2500' pcc4n5p2500' lrcc4n5p2500' dq1n5p2500' pdq1n5p2500' dq4n5p2500' pdq4n5p2500' tsvqn5p2500' pvqrn5p2500' BFucn5p2500' BFind1n5p2500' BFcc1n5p2500' BFind4n5p2500' BFcc4n5p2500' BFvqn5p2500' BFdq1n5p2500' BFdq14n5p2500'];% BFdq2n5p2500']; 
%Mp52500=[pucn5p2500' lrucn5p2500' pindn5p2500' lrindn5p2500' pccn5p2500' lrccn5p2500' pind4n5p2500' lrind4n5p2500' pcc4n5p2500' lrcc4n5p2500' dq1n5p2500' pdq1n5p2500' dq4n5p2500' pdq4n5p2500' tsvqn5p2500' pvqrn5p2500' BFucn5p2500' BFind1n5p2500' BFcc1n5p2500' BFind4n5p2500' BFcc4n5p2500' BFvqn5p2500']; 
csvwrite('HS2500p_5final3.csv',Mp52500)
Mp12500=[hits1p2500' pucn1p2500' lrucn1p2500' pindn1p2500' lrindn1p2500' pccn1p2500' lrccn1p2500' pind4n1p2500' lrind4n1p2500' pcc4n1p2500' lrcc4n1p2500' dq1n1p2500' pdq1n1p2500' dq4n1p2500' pdq4n1p2500' tsvqn1p2500' pvqrn1p2500' BFucn1p2500' BFind1n1p2500' BFcc1n1p2500' BFind4n1p2500' BFcc4n1p2500' BFvqn1p2500' BFdq1n1p2500' BFdq14n1p2500'];% BFdq2n1p2500']; 
%Mp12500=[pucn1p2500' lrucn1p2500' pindn1p2500' lrindn1p2500' pccn1p2500' lrccn1p2500' pind4n1p2500' lrind4n1p2500' pcc4n1p2500' lrcc4n1p2500' dq1n1p2500' pdq1n1p2500' dq4n1p2500' pdq4n1p2500' tsvqn1p2500' pvqrn1p2500' BFucn1p2500' BFind1n1p2500' BFcc1n1p2500' BFind4n1p2500' BFcc4n1p2500' BFvqn1p2500']; 
csvwrite('HS2500p_1final3.csv',Mp12500)

n=1000;
for j=1:20000
j
[y, sig, VaR5, VaR1, VaR5HS, VaR1HS]=gengarchnHS2(n,a0,a1,b1);
  hits5s1000(j)=sum(y<VaR5);hits5p1000(j)=sum(y<VaR5HS);
  hits1s1000(j)=sum(y<VaR1);hits1p1000(j)=sum(y<VaR1HS);
  [pucn1s1000(j),lrucn1s1000(j),pccn1s1000(j),lrccn1s1000(j),pindn1s1000(j),lrindn1s1000(j),pcc4n1s1000(j),lrcc4n1s1000(j),pind4n1s1000(j),lrind4n1s1000(j),pdq1n1s1000(j),dq1n1s1000(j),pdq4n1s1000(j),dq4n1s1000(j),tsvqn1s1000(j),pvqrn1s1000(j)]=runucccdqvqtests(y,VaR1,0.01);
  [pucn5s1000(j),lrucn5s1000(j),pccn5s1000(j),lrccn5s1000(j),pindn5s1000(j),lrindn5s1000(j),pcc4n5s1000(j),lrcc4n5s1000(j),pind4n5s1000(j),lrind4n5s1000(j),pdq1n5s1000(j),dq1n5s1000(j),pdq4n5s1000(j),dq4n5s1000(j),tsvqn5s1000(j),pvqrn5s1000(j)]=runucccdqvqtests(y,VaR5,0.05);
  [pucn1p1000(j),lrucn1p1000(j),pccn1p1000(j),lrccn1p1000(j),pindn1p1000(j),lrindn1p1000(j),pcc4n1p1000(j),lrcc4n1p1000(j),pind4n1p1000(j),lrind4n1p1000(j),pdq1n1p1000(j),dq1n1p1000(j),pdq4n1p1000(j),dq4n1p1000(j),tsvqn1p1000(j),pvqrn1p1000(j)]=runucccdqvqtests(y,VaR1HS,0.01);
  [pucn5p1000(j),lrucn5p1000(j),pccn5p1000(j),lrccn5p1000(j),pindn5p1000(j),lrindn5p1000(j),pcc4n5p1000(j),lrcc4n5p1000(j),pind4n5p1000(j),lrind4n5p1000(j),pdq1n5p1000(j),dq1n5p1000(j),pdq4n5p1000(j),dq4n5p1000(j),tsvqn5p1000(j),pvqrn5p1000(j)]=runucccdqvqtests(y,VaR5HS,0.05);
  [BFucn1s1000(j),BFcc1n1s1000(j),BFind1n1s1000(j),BFcc4n1s1000(j),BFind4n1s1000(j),BFvqn1s1000(j)]=runBFtests(y,VaR1,0.01);
  [BFucn5s1000(j),BFcc1n5s1000(j),BFind1n5s1000(j),BFcc4n5s1000(j),BFind4n5s1000(j),BFvqn5s1000(j)]=runBFtests(y,VaR5,0.05);
 [BFucn1p1000(j),BFcc1n1p1000(j),BFind1n1p1000(j),BFcc4n1p1000(j),BFind4n1p1000(j),BFvqn1p1000(j)]=runBFtests(y,VaR1HS,0.01);
 [BFucn5p1000(j),BFcc1n5p1000(j),BFind1n5p1000(j),BFcc4n5p1000(j),BFind4n5p1000(j),BFvqn5p1000(j)]=runBFtests(y,VaR5HS,0.05);
  BFdq1n1s1000(j) = 1/bfdqmlik3(y,VaR1,0.01);BFdq1n5s1000(j)=1/bfdqmlik3(y,VaR5,0.05);
  BFdq1n1p1000(j) = 1/bfdqmlik3(y,VaR1HS,0.01);BFdq1n5p1000(j)=1/bfdqmlik3(y,VaR5HS,0.05);
  BFdq14n1s1000(j) = 1/bfdqmlik4(y,VaR1,0.01);BFdq14n5s1000(j)=1/bfdqmlik4(y,VaR5,0.05);
  BFdq14n1p1000(j) = 1/bfdqmlik4(y,VaR1HS,0.01);BFdq14n5p1000(j)=1/bfdqmlik4(y,VaR5HS,0.05);
  BFdq2n1s1000(j) = 1/bfdq2mlik2(y,VaR1,0.01);BFdq2n5s1000(j)=1/bfdq2mlik2(y,VaR5,0.05);
  BFdq2n1p1000(j) = 1/bfdq2mlik2(y,VaR1HS,0.01);BFdq2n5p1000(j)=1/bfdq2mlik2(y,VaR5HS,0.05);
  
  if rem(j,100)==0
     j
     [pucn1s1000(j),lrucn1s1000(j),pccn1s1000(j),lrccn1s1000(j),pindn1s1000(j),lrindn1s1000(j),pcc4n1s1000(j),lrcc4n1s1000(j),pind4n1s1000(j),lrind4n1s1000(j),pdq1n1s1000(j),dq1n1s1000(j),pdq4n1s1000(j),dq4n1s1000(j),tsvqn1s1000(j),pvqrn1s1000(j);
      pucn1p1000(j),lrucn1p1000(j),pccn1p1000(j),lrccn1p1000(j),pindn1p1000(j),lrindn1p1000(j),pcc4n1p1000(j),lrcc4n1p1000(j),pind4n1p1000(j),lrind4n1p1000(j),pdq1n1p1000(j),dq1n1p1000(j),pdq4n1p1000(j),dq4n1p1000(j),tsvqn1p1000(j),pvqrn1p1000(j)]
     %[2*log(BFdq2n1s1000(j)) 2*log(BFdq2n5s1000(j)) 2*log(BFdq14n1s1000(j)) 2*log(BFdq14n5s1000(j));
     %  2*log(BFdq2n1p1000(j)) 2*log(BFdq2n5p1000(j)) 2*log(BFdq14n1p1000(j)) 2*log(BFdq14n5p1000(j))]  
  end   
  
end  

'size'
[length(pucn5s1000(pucn5s1000<0.05)) length(pindn5s1000(pindn5s1000<0.05)) length(pccn5s1000(pccn5s1000<0.05)) length(pind4n5s1000(pind4n5s1000<0.05)) length(pcc4n5s1000(pcc4n5s1000<0.05)) length(pdq1n5s1000(pdq1n5s1000<0.05)) length(pdq4n5s1000(pdq4n5s1000<0.05)) length(pvqrn5s1000(pvqrn5s1000<0.05)) length(BFucn5s1000(BFucn5s1000>1)) length(BFind1n5s1000(BFind1n5s1000>1)) length(BFcc1n5s1000(BFcc1n5s1000>1)) length(BFind4n5s1000(BFind4n5s1000>1)) length(BFcc4n5s1000(BFcc4n5s1000>1)) length(BFvqn5s1000(BFvqn5s1000>1)) length(BFdq1n5s1000(BFdq1n5s1000>1)) length(BFdq14n5s1000(BFdq14n5s1000>1)) length(BFdq2n5s1000(BFdq2n5s1000>1))]./length(pvqrn5s1000)
[length(pucn1s1000(pucn1s1000<0.05)) length(pindn1s1000(pindn1s1000<0.05)) length(pccn1s1000(pccn1s1000<0.05)) length(pind4n1s1000(pind4n1s1000<0.05)) length(pcc4n1s1000(pcc4n1s1000<0.05)) length(pdq1n1s1000(pdq1n1s1000<0.05)) length(pdq4n1s1000(pdq4n1s1000<0.05)) length(pvqrn1s1000(pvqrn1s1000<0.05)) length(BFucn1s1000(BFucn1s1000>1)) length(BFind1n1s1000(BFind1n1s1000>1)) length(BFcc1n1s1000(BFcc1n1s1000>1))  length(BFind4n1s1000(BFind4n1s1000>1)) length(BFcc4n1s1000(BFcc4n1s1000>1)) length(BFvqn1s1000(BFvqn1s1000>1)) length(BFdq1n1s1000(BFdq1n1s1000>1)) length(BFdq14n1s1000(BFdq14n1s1000>1)) length(BFdq2n1s1000(BFdq2n1s1000>1))]./length(pvqrn1s1000)
%[length(pucn5s1000(pucn5s1000<0.05)) length(pindn5s1000(pindn5s1000<0.05)) length(pccn5s1000(pccn5s1000<0.05)) length(pind4n5s1000(pind4n5s1000<0.05)) length(pcc4n5s1000(pcc4n5s1000<0.05)) length(pdq1n5s1000(pdq1n5s1000<0.05)) length(pdq4n5s1000(pdq4n5s1000<0.05)) length(pvqrn5s1000(pvqrn5s1000<0.05)) length(BFucn5s1000(BFucn5s1000>1)) length(BFind1n5s1000(BFind1n5s1000>1)) length(BFcc1n5s1000(BFcc1n5s1000>1)) length(BFind4n5s1000(BFind4n5s1000>1)) length(BFcc4n5s1000(BFcc4n5s1000>1)) length(BFvqn5s1000(BFvqn5s1000>1)) length(BFdq1n5s1000(BFdq1n5s1000>1)) length(BFdq14n5s1000(BFdq14n5s1000>1))]./length(pvqrn5s1000)
%[length(pucn1s1000(pucn1s1000<0.05)) length(pindn1s1000(pindn1s1000<0.05)) length(pccn1s1000(pccn1s1000<0.05)) length(pind4n1s1000(pind4n1s1000<0.05)) length(pcc4n1s1000(pcc4n1s1000<0.05)) length(pdq1n1s1000(pdq1n1s1000<0.05)) length(pdq4n1s1000(pdq4n1s1000<0.05)) length(pvqrn1s1000(pvqrn1s1000<0.05)) length(BFucn1s1000(BFucn1s1000>1)) length(BFind1n1s1000(BFind1n1s1000>1)) length(BFcc1n1s1000(BFcc1n1s1000>1))  length(BFind4n1s1000(BFind4n1s1000>1)) length(BFcc4n1s1000(BFcc4n1s1000>1)) length(BFvqn1s1000(BFvqn1s1000>1)) length(BFdq1n1s1000(BFdq1n1s1000>1)) length(BFdq14n1s1000(BFdq14n1s1000>1))]./length(pvqrn1s1000)


'power'
[length(pucn5p1000(pucn5p1000<0.05)) length(pindn5p1000(pindn5p1000<0.05)) length(pccn5p1000(pccn5p1000<0.05)) length(pind4n5p1000(pind4n5p1000<0.05)) length(pcc4n5p1000(pcc4n5p1000<0.05)) length(pdq1n5p1000(pdq1n5p1000<0.05)) length(pdq4n5p1000(pdq4n5p1000<0.05)) length(pvqrn5p1000(pvqrn5p1000<0.05)) length(BFucn5p1000(BFucn5p1000>1)) length(BFind1n5p1000(BFind1n5p1000>1)) length(BFcc1n5p1000(BFcc1n5p1000>1)) length(BFind4n5p1000(BFind4n5p1000>1)) length(BFcc4n5p1000(BFcc4n5p1000>1)) length(BFvqn5p1000(BFvqn5p1000>1)) length(BFdq1n5p1000(BFdq1n5p1000>1)) length(BFdq14n5p1000(BFdq14n5p1000>1)) length(BFdq2n5p1000(BFdq2n5p1000>1))]./length(pvqrn5p1000)
[length(pucn1p1000(pucn1p1000<0.05)) length(pindn1p1000(pindn1p1000<0.05)) length(pccn1p1000(pccn1p1000<0.05)) length(pind4n1p1000(pind4n1p1000<0.05)) length(pcc4n1p1000(pcc4n1p1000<0.05)) length(pdq1n1p1000(pdq1n1p1000<0.05)) length(pdq4n1p1000(pdq4n1p1000<0.05)) length(pvqrn1p1000(pvqrn1p1000<0.05)) length(BFucn1p1000(BFucn1p1000>1)) length(BFind1n1p1000(BFind1n1p1000>1)) length(BFcc1n1p1000(BFcc1n1p1000>1)) length(BFind4n1p1000(BFind4n1p1000>1)) length(BFcc4n1p1000(BFcc4n1p1000>1)) length(BFvqn1p1000(BFvqn1p1000>1)) length(BFdq1n1p1000(BFdq1n1p1000>1)) length(BFdq14n1p1000(BFdq14n5p1000>1)) length(BFdq2n1p1000(BFdq2n1p1000>1))]./length(pvqrn1p1000)
%[length(pucn5p1000(pucn5p1000<0.05)) length(pindn5p1000(pindn5p1000<0.05)) length(pccn5p1000(pccn5p1000<0.05)) length(pind4n5p1000(pind4n5p1000<0.05)) length(pcc4n5p1000(pcc4n5p1000<0.05)) length(pdq1n5p1000(pdq1n5p1000<0.05)) length(pdq4n5p1000(pdq4n5p1000<0.05)) length(pvqrn5p1000(pvqrn5p1000<0.05)) length(BFucn5p1000(BFucn5p1000>1)) length(BFind1n5p1000(BFind1n5p1000>1)) length(BFcc1n5p1000(BFcc1n5p1000>1)) length(BFind4n5p1000(BFind4n5p1000>1)) length(BFcc4n5p1000(BFcc4n5p1000>1)) length(BFvqn5p1000(BFvqn5p1000>1)) length(BFdq1n5p1000(BFdq1n5p1000>1)) length(BFdq14n5p1000(BFdq14n5p1000>1))]./length(pvqrn5p1000)
%[length(pucn1p1000(pucn1p1000<0.05)) length(pindn1p1000(pindn1p1000<0.05)) length(pccn1p1000(pccn1p1000<0.05)) length(pind4n1p1000(pind4n1p1000<0.05)) length(pcc4n1p1000(pcc4n1p1000<0.05)) length(pdq1n1p1000(pdq1n1p1000<0.05)) length(pdq4n1p1000(pdq4n1p1000<0.05)) length(pvqrn1p1000(pvqrn1p1000<0.05)) length(BFucn1p1000(BFucn1p1000>1)) length(BFind1n1p1000(BFind1n1p1000>1)) length(BFcc1n1p1000(BFcc1n1p1000>1)) length(BFind4n1p1000(BFind4n1p1000>1)) length(BFcc4n1p1000(BFcc4n1p1000>1)) length(BFvqn1p1000(BFvqn1p1000>1)) length(BFdq1n1p1000(BFdq1n1p1000>1)) length(BFdq14n1p1000(BFdq14n1p1000>1))]./length(pvqrn1p1000)

BFuc5s1000c=prctile(BFucn5s1000,95);BFind5s1000c=prctile(BFind1n5s1000,95);
BFcc5s1000c=prctile(BFcc1n5s1000,95);BFcc45s1000c=prctile(BFcc4n5s1000,95);
BFind45s1000c=prctile(BFind4n5s1000,95);BFvq5s1000c=prctile(BFvqn5s1000,95);
BFdq15s1000c=prctile(BFdq1n5s1000,95);BFdq145s1000c=prctile(BFdq14n5s1000,95);BFdq25s1000c=prctile(BFdq2n5s1000,95);
uc5s1000c=prctile(pucn5s1000,5);ind5s1000c=prctile(pindn5s1000,5);ind45s1000c=prctile(pind4n5s1000,5);
cc5s1000c=prctile(pccn5s1000,5);cc45s1000c=prctile(pcc4n5s1000,5);
dq15s1000c=prctile(dq1n5s1000,95);dq45s1000c=prctile(dq4n5s1000,95);tsvq5s1000c=prctile(tsvqn5s1000,95);

BFuc1s1000c=prctile(BFucn1s1000,95);BFind1s1000c=prctile(BFind1n1s1000,95);
BFcc1s1000c=prctile(BFcc1n1s1000,95);BFcc41s1000c=prctile(BFcc4n1s1000,95);
BFind41s1000c=prctile(BFind4n1s1000,95);BFvq1s1000c=prctile(BFvqn1s1000,95);
BFdq11s1000c=prctile(BFdq1n1s1000,95);BFdq141s1000c=prctile(BFdq14n1s1000,95);BFdq21s1000c=prctile(BFdq2n1s1000,95);
uc1s1000c=prctile(pucn1s1000,5);ind1s1000c=prctile(pindn1s1000,5);ind41s1000c=prctile(pind4n1s1000,5);
cc1s1000c=prctile(pccn1s1000,5);cc41s1000c=prctile(pcc4n1s1000,5);
dq11s1000c=prctile(dq1n1s1000,95);dq41s1000c=prctile(dq4n1s1000,95);tsvq1s1000c=prctile(tsvqn1s1000,95);

disp('corrected size');
[length(pucn5s1000(pucn5s1000<=uc5s1000c)) length(pindn5s1000(pindn5s1000<ind5s1000c)) length(pccn5s1000(pccn5s1000<cc5s1000c)) length(pind4n5s1000(pind4n5s1000<ind45s1000c)) length(pcc4n5s1000(pcc4n5s1000<cc45s1000c)) length(dq1n5s1000(dq1n5s1000>dq15s1000c)) length(dq4n5s1000(dq4n5s1000>dq45s1000c)) length(tsvqn5s1000(tsvqn5s1000>tsvq5s1000c)) length(BFucn5s1000(BFucn5s1000>BFuc5s1000c)) length(BFind1n5s1000(BFind1n5s1000>BFind5s1000c)) length(BFcc1n5s1000(BFcc1n5s1000>BFcc5s1000c)) length(BFind4n5s1000(BFind4n5s1000>BFind45s1000c)) length(BFcc4n5s1000(BFcc4n5s1000>BFcc45s1000c)) length(BFvqn5s1000(BFvqn5s1000>BFvq5s1000c)) length(BFdq1n5s1000(BFdq1n5s1000>BFdq15s1000c)) length(BFdq14n5s1000(BFdq14n5s1000>BFdq145s1000c)) length(BFdq2n5s1000(BFdq2n5s1000>BFdq25s1000c))]./length(pvqrn5s1000)
[length(pucn1s1000(pucn1s1000<=uc1s1000c)) length(pindn1s1000(pindn1s1000<ind1s1000c)) length(pccn1s1000(pccn1s1000<cc1s1000c)) length(pind4n1s1000(pind4n1s1000<ind41s1000c)) length(pcc4n1s1000(pcc4n1s1000<cc41s1000c)) length(dq1n1s1000(dq1n1s1000>dq11s1000c)) length(dq4n1s1000(dq4n1s1000>dq41s1000c)) length(tsvqn1s1000(tsvqn1s1000>tsvq1s1000c)) length(BFucn1s1000(BFucn1s1000>BFuc1s1000c)) length(BFind1n1s1000(BFind1n1s1000>BFind1s1000c)) length(BFcc1n1s1000(BFcc1n1s1000>BFcc1s1000c)) length(BFind4n1s1000(BFind4n1s1000>BFind41s1000c)) length(BFcc4n1s1000(BFcc4n1s1000>BFcc41s1000c)) length(BFvqn1s1000(BFvqn1s1000>BFvq1s1000c)) length(BFdq1n1s1000(BFdq1n1s1000>BFdq11s1000c)) length(BFdq14n1s1000(BFdq14n1s1000>BFdq141s1000c)) length(BFdq2n1s1000(BFdq2n1s1000>BFdq21s1000c))]./length(pvqrn1s1000)
%[length(pucn5s1000(pucn5s1000<=uc5s1000c)) length(pindn5s1000(pindn5s1000<ind5s1000c)) length(pccn5s1000(pccn5s1000<cc5s1000c)) length(pind4n5s1000(pind4n5s1000<ind45s1000c)) length(pcc4n5s1000(pcc4n5s1000<cc45s1000c)) length(dq1n5s1000(dq1n5s1000>dq15s1000c)) length(dq4n5s1000(dq4n5s1000>dq45s1000c)) length(tsvqn5s1000(tsvqn5s1000>tsvq5s1000c)) length(BFucn5s1000(BFucn5s1000>BFuc5s1000c)) length(BFind1n5s1000(BFind1n5s1000>BFind5s1000c)) length(BFcc1n5s1000(BFcc1n5s1000>BFcc5s1000c)) length(BFind4n5s1000(BFind4n5s1000>BFind45s1000c)) length(BFcc4n5s1000(BFcc4n5s1000>BFcc45s1000c)) length(BFvqn5s1000(BFvqn5s1000>BFvq5s1000c)) length(BFdq1n5s1000(BFdq1n5s1000>BFdq15s1000c)) length(BFdq14n5s1000(BFdq14n5s1000>BFdq145s1000c))]./length(pvqrn5s1000)
%[length(pucn1s1000(pucn1s1000<=uc1s1000c)) length(pindn1s1000(pindn1s1000<ind1s1000c)) length(pccn1s1000(pccn1s1000<cc1s1000c)) length(pind4n1s1000(pind4n1s1000<ind41s1000c)) length(pcc4n1s1000(pcc4n1s1000<cc41s1000c)) length(dq1n1s1000(dq1n1s1000>dq11s1000c)) length(dq4n1s1000(dq4n1s1000>dq41s1000c)) length(tsvqn1s1000(tsvqn1s1000>tsvq1s1000c)) length(BFucn1s1000(BFucn1s1000>BFuc1s1000c)) length(BFind1n1s1000(BFind1n1s1000>BFind1s1000c)) length(BFcc1n1s1000(BFcc1n1s1000>BFcc1s1000c)) length(BFind4n1s1000(BFind4n1s1000>BFind41s1000c)) length(BFcc4n1s1000(BFcc4n1s1000>BFcc41s1000c)) length(BFvqn1s1000(BFvqn1s1000>BFvq1s1000c)) length(BFdq1n1s1000(BFdq1n1s1000>BFdq11s1000c)) length(BFdq14n1s1000(BFdq14n1s1000>BFdq141s1000c))]./length(pvqrn1s1000)


disp('size-adjusted power');
[length(pucn5p1000(pucn5p1000<=uc5s1000c)) length(pindn5p1000(pindn5p1000<ind5s1000c)) length(pccn5p1000(pccn5p1000<cc5s1000c)) length(pind4n5p1000(pind4n5p1000<ind45s1000c)) length(pcc4n5p1000(pcc4n5p1000<cc45s1000c)) length(dq1n5p1000(dq1n5p1000>dq15s1000c)) length(dq4n5p1000(dq4n5p1000>dq45s1000c)) length(tsvqn5p1000(tsvqn5p1000>tsvq5s1000c)) length(BFucn5p1000(BFucn5p1000>BFuc5s1000c)) length(BFind1n5p1000(BFind1n5p1000>BFind5s1000c)) length(BFcc1n5p1000(BFcc1n5p1000>BFcc5s1000c)) length(BFind4n5p1000(BFind4n5p1000>BFind45s1000c)) length(BFcc4n5p1000(BFcc4n5p1000>BFcc45s1000c)) length(BFvqn5p1000(BFvqn5p1000>BFvq5s1000c)) length(BFdq1n5p1000(BFdq1n5p1000>BFdq15s1000c)) length(BFdq14n5p1000(BFdq14n5p1000>BFdq145s1000c)) length(BFdq2n5p1000(BFdq2n5p1000>BFdq25s1000c))]./length(pvqrn5p1000)
[length(pucn1p1000(pucn1p1000<=uc1s1000c)) length(pindn1p1000(pindn1p1000<ind1s1000c)) length(pccn1p1000(pccn1p1000<cc1s1000c)) length(pind4n1p1000(pind4n1p1000<ind41s1000c)) length(pcc4n1p1000(pcc4n1p1000<cc41s1000c)) length(dq1n1p1000(dq1n1p1000>dq11s1000c)) length(dq4n1p1000(dq4n1p1000>dq41s1000c)) length(tsvqn1p1000(tsvqn1p1000>tsvq1s1000c)) length(BFucn1p1000(BFucn1p1000>BFuc1s1000c)) length(BFind1n1p1000(BFind1n1p1000>BFind1s1000c)) length(BFcc1n1p1000(BFcc1n1p1000>BFcc1s1000c)) length(BFind4n1p1000(BFind4n1p1000>BFind41s1000c)) length(BFcc4n1p1000(BFcc4n1p1000>BFcc41s1000c)) length(BFvqn1p1000(BFvqn1p1000>BFvq1s1000c)) length(BFdq1n1p1000(BFdq1n1p1000>BFdq11s1000c)) length(BFdq14n1p1000(BFdq14n1p1000>BFdq141s1000c)) length(BFdq2n1p1000(BFdq2n1p1000>BFdq21s1000c))]./length(pvqrn1p1000)
%[length(pucn5p1000(pucn5p1000<=uc5s1000c)) length(pindn5p1000(pindn5p1000<ind5s1000c)) length(pccn5p1000(pccn5p1000<cc5s1000c)) length(pind4n5p1000(pind4n5p1000<ind45s1000c)) length(pcc4n5p1000(pcc4n5p1000<cc45s1000c)) length(dq1n5p1000(dq1n5p1000>dq15s1000c)) length(dq4n5p1000(dq4n5p1000>dq45s1000c)) length(tsvqn5p1000(tsvqn5p1000>tsvq5s1000c)) length(BFucn5p1000(BFucn5p1000>BFuc5s1000c)) length(BFind1n5p1000(BFind1n5p1000>BFind5s1000c)) length(BFcc1n5p1000(BFcc1n5p1000>BFcc5s1000c)) length(BFind4n5p1000(BFind4n5p1000>BFind45s1000c)) length(BFcc4n5p1000(BFcc4n5p1000>BFcc45s1000c)) length(BFvqn5p1000(BFvqn5p1000>BFvq5s1000c)) length(BFdq1n5p1000(BFdq1n5p1000>BFdq15s1000c)) length(BFdq14n5p1000(BFdq14n5p1000>BFdq145s1000c))]./length(pvqrn5p1000)
%[length(pucn1p1000(pucn1p1000<=uc1s1000c)) length(pindn1p1000(pindn1p1000<ind1s1000c)) length(pccn1p1000(pccn1p1000<cc1s1000c)) length(pind4n1p1000(pind4n1p1000<ind41s1000c)) length(pcc4n1p1000(pcc4n1p1000<cc41s1000c)) length(dq1n1p1000(dq1n1p1000>dq11s1000c)) length(dq4n1p1000(dq4n1p1000>dq41s1000c)) length(tsvqn1p1000(tsvqn1p1000>tsvq1s1000c)) length(BFucn1p1000(BFucn1p1000>BFuc1s1000c)) length(BFind1n1p1000(BFind1n1p1000>BFind1s1000c)) length(BFcc1n1p1000(BFcc1n1p1000>BFcc1s1000c)) length(BFind4n1p1000(BFind4n1p1000>BFind41s1000c)) length(BFcc4n1p1000(BFcc4n1p1000>BFcc41s1000c)) length(BFvqn1p1000(BFvqn1p1000>BFvq1s1000c)) length(BFdq1n1p1000(BFdq1n1p1000>BFdq11s1000c)) length(BFdq14n1p1000(BFdq14n1p1000>BFdq141s1000c))]./length(pvqrn1p1000)

Ms51000=[hits5s1000' pucn5s1000' lrucn5s1000' pindn5s1000' lrindn5s1000' pccn5s1000' lrccn5s1000' pind4n5s1000' lrind4n5s1000' pcc4n5s1000' lrcc4n5s1000' dq1n5s1000' pdq1n5s1000' dq4n5s1000' pdq4n5s1000' tsvqn5s1000' pvqrn5s1000' BFucn5s1000' BFind1n5s1000' BFcc1n5s1000' BFind4n5s1000' BFcc4n5s1000' BFvqn5s1000' BFdq1n5s1000' BFdq14n5s1000' BFdq2n5s1000']; 
%Ms51000=[hits5s1000' pucn5s1000' lrucn5s1000' pindn5s1000' lrindn5s1000' pccn5s1000' lrccn5s1000' pind4n5s1000' lrind4n5s1000' pcc4n5s1000' lrcc4n5s1000' dq1n5s1000' pdq1n5s1000' dq4n5s1000' pdq4n5s1000' tsvqn5s1000' pvqrn5s1000' BFucn5s1000' BFind1n5s1000' BFcc1n5s1000' BFind4n5s1000' BFcc4n5s1000' BFvqn5s1000']; 
csvwrite('HS1000s_5final3.csv',Ms51000)
Ms11000=[hits1s1000' pucn1s1000' lrucn1s1000' pindn1s1000' lrindn1s1000' pccn1s1000' lrccn1s1000' pind4n1s1000' lrind4n1s1000' pcc4n1s1000' lrcc4n1s1000' dq1n1s1000' pdq1n1s1000' dq4n1s1000' pdq4n1s1000' tsvqn1s1000' pvqrn1s1000' BFucn1s1000' BFind1n1s1000' BFcc1n1s1000' BFind4n1s1000' BFcc4n1s1000' BFvqn1s1000' BFdq1n1s1000' BFdq14n1s1000' BFdq2n1s1000']; 
%Ms11000=[hits1s1000' pucn1s1000' lrucn1s1000' pindn1s1000' lrindn1s1000' pccn1s1000' lrccn1s1000' pind4n1s1000' lrind4n1s1000' pcc4n1s1000' lrcc4n1s1000' dq1n1s1000' pdq1n1s1000' dq4n1s1000' pdq4n1s1000' tsvqn1s1000' pvqrn1s1000' BFucn1s1000' BFind1n1s1000' BFcc1n1s1000' BFind4n1s1000' BFcc4n1s1000' BFvqn1s1000']; 
csvwrite('HS1000s_1final3.csv',Ms11000)
Mp51000=[hits5p1000' pucn5p1000' lrucn5p1000' pindn5p1000' lrindn5p1000' pccn5p1000' lrccn5p1000' pind4n5p1000' lrind4n5p1000' pcc4n5p1000' lrcc4n5p1000' dq1n5p1000' pdq1n5p1000' dq4n5p1000' pdq4n5p1000' tsvqn5p1000' pvqrn5p1000' BFucn5p1000' BFind1n5p1000' BFcc1n5p1000' BFind4n5p1000' BFcc4n5p1000' BFvqn5p1000' BFdq1n5p1000' BFdq14n5p1000' BFdq2n5p1000']; 
%Mp51000=[hits5p1000' pucn5p1000' lrucn5p1000' pindn5p1000' lrindn5p1000' pccn5p1000' lrccn5p1000' pind4n5p1000' lrind4n5p1000' pcc4n5p1000' lrcc4n5p1000' dq1n5p1000' pdq1n5p1000' dq4n5p1000' pdq4n5p1000' tsvqn5p1000' pvqrn5p1000' BFucn5p1000' BFind1n5p1000' BFcc1n5p1000' BFind4n5p1000' BFcc4n5p1000' BFvqn5p1000']; 
csvwrite('HS1000p_5final3.csv',Mp51000)
Mp11000=[hits1p1000' pucn1p1000' lrucn1p1000' pindn1p1000' lrindn1p1000' pccn1p1000' lrccn1p1000' pind4n1p1000' lrind4n1p1000' pcc4n1p1000' lrcc4n1p1000' dq1n1p1000' pdq1n1p1000' dq4n1p1000' pdq4n1p1000' tsvqn1p1000' pvqrn1p1000' BFucn1p1000' BFind1n1p1000' BFcc1n1p1000' BFind4n1p1000' BFcc4n1p1000' BFvqn1p1000' BFdq1n1p1000' BFdq14n1p1000' BFdq2n1p1000']; 
%Mp11000=[hits1p1000' pucn1p1000' lrucn1p1000' pindn1p1000' lrindn1p1000' pccn1p1000' lrccn1p1000' pind4n1p1000' lrind4n1p1000' pcc4n1p1000' lrcc4n1p1000' dq1n1p1000' pdq1n1p1000' dq4n1p1000' pdq4n1p1000' tsvqn1p1000' pvqrn1p1000' BFucn1p1000' BFind1n1p1000' BFcc1n1p1000' BFind4n1p1000' BFcc4n1p1000' BFvqn1p1000']; 
csvwrite('HS1000p_1final3.csv',Mp11000)

n=500;
for j=1:5000
j
[y, sig, VaR5, VaR1, VaR5HS, VaR1HS]=gengarchnHS2(n,a0,a1,b1);
  hits5s500(j)=sum(y<VaR5);hits5p500(j)=sum(y<VaR5HS);
  hits1s500(j)=sum(y<VaR1);hits1p500(j)=sum(y<VaR1HS);
  [pucn1s500(j),lrucn1s500(j),pccn1s500(j),lrccn1s500(j),pindn1s500(j),lrindn1s500(j),pcc4n1s500(j),lrcc4n1s500(j),pind4n1s500(j),lrind4n1s500(j),pdq1n1s500(j),dq1n1s500(j),pdq4n1s500(j),dq4n1s500(j),tsvqn1s500(j),pvqrn1s500(j)]=runucccdqvqtests(y,VaR1,0.01);
  [pucn5s500(j),lrucn5s500(j),pccn5s500(j),lrccn5s500(j),pindn5s500(j),lrindn5s500(j),pcc4n5s500(j),lrcc4n5s500(j),pind4n5s500(j),lrind4n5s500(j),pdq1n5s500(j),dq1n5s500(j),pdq4n5s500(j),dq4n5s500(j),tsvqn5s500(j),pvqrn5s500(j)]=runucccdqvqtests(y,VaR5,0.05);
  [pucn1p500(j),lrucn1p500(j),pccn1p500(j),lrccn1p500(j),pindn1p500(j),lrindn1p500(j),pcc4n1p500(j),lrcc4n1p500(j),pind4n1p500(j),lrind4n1p500(j),pdq1n1p500(j),dq1n1p500(j),pdq4n1p500(j),dq4n1p500(j),tsvqn1p500(j),pvqrn1p500(j)]=runucccdqvqtests(y,VaR1HS,0.01);
  [pucn5p500(j),lrucn5p500(j),pccn5p500(j),lrccn5p500(j),pindn5p500(j),lrindn5p500(j),pcc4n5p500(j),lrcc4n5p500(j),pind4n5p500(j),lrind4n5p500(j),pdq1n5p500(j),dq1n5p500(j),pdq4n5p500(j),dq4n5p500(j),tsvqn5p500(j),pvqrn5p500(j)]=runucccdqvqtests(y,VaR5HS,0.05);
  [BFucn1s500(j),BFcc1n1s500(j),BFind1n1s500(j),BFcc4n1s500(j),BFind4n1s500(j),BFvqn1s500(j)]=runBFtests(y,VaR1,0.01);
  [BFucn5s500(j),BFcc1n5s500(j),BFind1n5s500(j),BFcc4n5s500(j),BFind4n5s500(j),BFvqn5s500(j)]=runBFtests(y,VaR5,0.05);
 [BFucn1p500(j),BFcc1n1p500(j),BFind1n1p500(j),BFcc4n1p500(j),BFind4n1p500(j),BFvqn1p500(j)]=runBFtests(y,VaR1HS,0.01);
 [BFucn5p500(j),BFcc1n5p500(j),BFind1n5p500(j),BFcc4n5p500(j),BFind4n5p500(j),BFvqn5p500(j)]=runBFtests(y,VaR5HS,0.05);
  BFdq1n1s500(j) = 1/bfdqmlik3(y,VaR1,0.01);BFdq1n5s500(j)=1/bfdqmlik3(y,VaR5,0.05);
  BFdq1n1p500(j) = 1/bfdqmlik3(y,VaR1HS,0.01);BFdq1n5p500(j)=1/bfdqmlik3(y,VaR5HS,0.05);
  BFdq14n1s500(j) = 1/bfdqmlik4(y,VaR1,0.01);BFdq14n5s500(j)=1/bfdqmlik4(y,VaR5,0.05);
  BFdq14n1p500(j) = 1/bfdqmlik4(y,VaR1HS,0.01);BFdq14n5p500(j)=1/bfdqmlik4(y,VaR5HS,0.05);
%   BFdq2n1s500(j) = 1/bfdq2mlik2(y,VaR1,0.01);BFdq2n5s500(j)=1/bfdq2mlik2(y,VaR5,0.05);
%   BFdq2n1p500(j) = 1/bfdq2mlik2(y,VaR1HS,0.01);BFdq2n5p500(j)=1/bfdq2mlik2(y,VaR5HS,0.05);
  
  if rem(j,100)==0
     j
     [pucn1s500(j),lrucn1s500(j),pccn1s500(j),lrccn1s500(j),pindn1s500(j),lrindn1s500(j),pcc4n1s500(j),lrcc4n1s500(j),pind4n1s500(j),lrind4n1s500(j),pdq1n1s500(j),dq1n1s500(j),pdq4n1s500(j),dq4n1s500(j),tsvqn1s500(j),pvqrn1s500(j);
      pucn1p500(j),lrucn1p500(j),pccn1p500(j),lrccn1p500(j),pindn1p500(j),lrindn1p500(j),pcc4n1p500(j),lrcc4n1p500(j),pind4n1p500(j),lrind4n1p500(j),pdq1n1p500(j),dq1n1p500(j),pdq4n1p500(j),dq4n1p500(j),tsvqn1p500(j),pvqrn1p500(j)]
     %[2*log(BFdq2n1s500(j)) 2*log(BFdq2n5s500(j)) 2*log(BFdq14n1s500(j)) 2*log(BFdq14n5s500(j));
     %  2*log(BFdq2n1p500(j)) 2*log(BFdq2n5p500(j)) 2*log(BFdq14n1p500(j)) 2*log(BFdq14n5p500(j))]  
  end   
  
end  

disp('size');
[length(pucn5s500(pucn5s500<0.05)) length(pindn5s500(pindn5s500<0.05)) length(pccn5s500(pccn5s500<0.05)) length(pind4n5s500(pind4n5s500<0.05)) length(pcc4n5s500(pcc4n5s500<0.05)) length(pdq1n5s500(pdq1n5s500<0.05)) length(pdq4n5s500(pdq4n5s500<0.05)) length(pvqrn5s500(pvqrn5s500<0.05)) length(BFucn5s500(BFucn5s500>1)) length(BFind1n5s500(BFind1n5s500>1)) length(BFcc1n5s500(BFcc1n5s500>1)) length(BFind4n5s500(BFind4n5s500>1)) length(BFcc4n5s500(BFcc4n5s500>1)) length(BFvqn5s500(BFvqn5s500>1)) length(BFdq1n5s500(BFdq1n5s500>1)) length(BFdq14n5s500(BFdq14n5s500>1)) length(BFdq2n5s500(BFdq2n5s500>1))]./length(pvqrn5s500)
[length(pucn1s500(pucn1s500<0.05)) length(pindn1s500(pindn1s500<0.05)) length(pccn1s500(pccn1s500<0.05)) length(pind4n1s500(pind4n1s500<0.05)) length(pcc4n1s500(pcc4n1s500<0.05)) length(pdq1n1s500(pdq1n1s500<0.05)) length(pdq4n1s500(pdq4n1s500<0.05)) length(pvqrn1s500(pvqrn1s500<0.05)) length(BFucn1s500(BFucn1s500>1)) length(BFind1n1s500(BFind1n1s500>1)) length(BFcc1n1s500(BFcc1n1s500>1))  length(BFind4n1s500(BFind4n1s500>1)) length(BFcc4n1s500(BFcc4n1s500>1)) length(BFvqn1s500(BFvqn1s500>1)) length(BFdq1n1s500(BFdq1n1s500>1)) length(BFdq14n1s500(BFdq14n1s500>1)) length(BFdq2n1s500(BFdq2n1s500>1))]./length(pvqrn1s500)
%[length(pucn5s500(pucn5s500<0.05)) length(pindn5s500(pindn5s500<0.05)) length(pccn5s500(pccn5s500<0.05)) length(pind4n5s500(pind4n5s500<0.05)) length(pcc4n5s500(pcc4n5s500<0.05)) length(pdq1n5s500(pdq1n5s500<0.05)) length(pdq4n5s500(pdq4n5s500<0.05)) length(pvqrn5s500(pvqrn5s500<0.05)) length(BFucn5s500(BFucn5s500>1)) length(BFind1n5s500(BFind1n5s500>1)) length(BFcc1n5s500(BFcc1n5s500>1)) length(BFind4n5s500(BFind4n5s500>1)) length(BFcc4n5s500(BFcc4n5s500>1)) length(BFvqn5s500(BFvqn5s500>1)) length(BFdq1n5s500(BFdq1n5s500>1)) length(BFdq14n5s500(BFdq14n5s500>1))]./length(pvqrn5s500)
%[length(pucn1s500(pucn1s500<0.05)) length(pindn1s500(pindn1s500<0.05)) length(pccn1s500(pccn1s500<0.05)) length(pind4n1s500(pind4n1s500<0.05)) length(pcc4n1s500(pcc4n1s500<0.05)) length(pdq1n1s500(pdq1n1s500<0.05)) length(pdq4n1s500(pdq4n1s500<0.05)) length(pvqrn1s500(pvqrn1s500<0.05)) length(BFucn1s500(BFucn1s500>1)) length(BFind1n1s500(BFind1n1s500>1)) length(BFcc1n1s500(BFcc1n1s500>1))  length(BFind4n1s500(BFind4n1s500>1)) length(BFcc4n1s500(BFcc4n1s500>1)) length(BFvqn1s500(BFvqn1s500>1)) length(BFdq1n1s500(BFdq1n1s500>1)) length(BFdq14n1s500(BFdq14n1s500>1))]./length(pvqrn1s500)


disp('power');
[length(pucn5p500(pucn5p500<0.05)) length(pindn5p500(pindn5p500<0.05)) length(pccn5p500(pccn5p500<0.05)) length(pind4n5p500(pind4n5p500<0.05)) length(pcc4n5p500(pcc4n5p500<0.05)) length(pdq1n5p500(pdq1n5p500<0.05)) length(pdq4n5p500(pdq4n5p500<0.05)) length(pvqrn5p500(pvqrn5p500<0.05)) length(BFucn5p500(BFucn5p500>1)) length(BFind1n5p500(BFind1n5p500>1)) length(BFcc1n5p500(BFcc1n5p500>1)) length(BFind4n5p500(BFind4n5p500>1)) length(BFcc4n5p500(BFcc4n5p500>1)) length(BFvqn5p500(BFvqn5p500>1)) length(BFdq1n5p500(BFdq1n5p500>1)) length(BFdq14n5p500(BFdq14n5p500>1)) length(BFdq2n5p500(BFdq2n5p500>1))]./length(pvqrn5p500)
[length(pucn1p500(pucn1p500<0.05)) length(pindn1p500(pindn1p500<0.05)) length(pccn1p500(pccn1p500<0.05)) length(pind4n1p500(pind4n1p500<0.05)) length(pcc4n1p500(pcc4n1p500<0.05)) length(pdq1n1p500(pdq1n1p500<0.05)) length(pdq4n1p500(pdq4n1p500<0.05)) length(pvqrn1p500(pvqrn1p500<0.05)) length(BFucn1p500(BFucn1p500>1)) length(BFind1n1p500(BFind1n1p500>1)) length(BFcc1n1p500(BFcc1n1p500>1)) length(BFind4n1p500(BFind4n1p500>1)) length(BFcc4n1p500(BFcc4n1p500>1)) length(BFvqn1p500(BFvqn1p500>1)) length(BFdq1n1p500(BFdq1n1p500>1)) length(BFdq14n1p500(BFdq14n5p500>1)) length(BFdq2n1p500(BFdq2n1p500>1))]./length(pvqrn1p500)
%[length(pucn5p500(pucn5p500<0.05)) length(pindn5p500(pindn5p500<0.05)) length(pccn5p500(pccn5p500<0.05)) length(pind4n5p500(pind4n5p500<0.05)) length(pcc4n5p500(pcc4n5p500<0.05)) length(pdq1n5p500(pdq1n5p500<0.05)) length(pdq4n5p500(pdq4n5p500<0.05)) length(pvqrn5p500(pvqrn5p500<0.05)) length(BFucn5p500(BFucn5p500>1)) length(BFind1n5p500(BFind1n5p500>1)) length(BFcc1n5p500(BFcc1n5p500>1)) length(BFind4n5p500(BFind4n5p500>1)) length(BFcc4n5p500(BFcc4n5p500>1)) length(BFvqn5p500(BFvqn5p500>1)) length(BFdq1n5p500(BFdq1n5p500>1)) length(BFdq14n5p500(BFdq14n5p500>1))]./length(pvqrn5p500)
%[length(pucn1p500(pucn1p500<0.05)) length(pindn1p500(pindn1p500<0.05)) length(pccn1p500(pccn1p500<0.05)) length(pind4n1p500(pind4n1p500<0.05)) length(pcc4n1p500(pcc4n1p500<0.05)) length(pdq1n1p500(pdq1n1p500<0.05)) length(pdq4n1p500(pdq4n1p500<0.05)) length(pvqrn1p500(pvqrn1p500<0.05)) length(BFucn1p500(BFucn1p500>1)) length(BFind1n1p500(BFind1n1p500>1)) length(BFcc1n1p500(BFcc1n1p500>1)) length(BFind4n1p500(BFind4n1p500>1)) length(BFcc4n1p500(BFcc4n1p500>1)) length(BFvqn1p500(BFvqn1p500>1)) length(BFdq1n1p500(BFdq1n1p500>1)) length(BFdq14n1p500(BFdq14n1p500>1))]./length(pvqrn1p500)

BFuc5s500c=prctile(BFucn5s500,95);BFind5s500c=prctile(BFind1n5s500,95);
BFcc5s500c=prctile(BFcc1n5s500,95);BFcc45s500c=prctile(BFcc4n5s500,95);
BFind45s500c=prctile(BFind4n5s500,95);BFvq5s500c=prctile(BFvqn5s500,95);
BFdq15s500c=prctile(BFdq1n5s500,95);BFdq145s500c=prctile(BFdq14n5s500,95);BFdq25s500c=prctile(BFdq2n5s500,95);
uc5s500c=prctile(pucn5s500,5);ind5s500c=prctile(pindn5s500,5);ind45s500c=prctile(pind4n5s500,5);
cc5s500c=prctile(pccn5s500,5);cc45s500c=prctile(pcc4n5s500,5);
dq15s500c=prctile(dq1n5s500,95);dq45s500c=prctile(dq4n5s500,95);tsvq5s500c=prctile(tsvqn5s500,95);

BFuc1s500c=prctile(BFucn1s500,95);BFind1s500c=prctile(BFind1n1s500,95);
BFcc1s500c=prctile(BFcc1n1s500,95);BFcc41s500c=prctile(BFcc4n1s500,95);
BFind41s500c=prctile(BFind4n1s500,95);BFvq1s500c=prctile(BFvqn1s500,95);
BFdq11s500c=prctile(BFdq1n1s500,95);BFdq141s500c=prctile(BFdq14n1s500,95);BFdq21s500c=prctile(BFdq2n1s500,95);
uc1s500c=prctile(pucn1s500,5);ind1s500c=prctile(pindn1s500,5);ind41s500c=prctile(pind4n1s500,5);
cc1s500c=prctile(pccn1s500,5);cc41s500c=prctile(pcc4n1s500,5);
dq11s500c=prctile(dq1n1s500,95);dq41s500c=prctile(dq4n1s500,95);tsvq1s500c=prctile(tsvqn1s500,95);

disp('corrected size');
[length(pucn5s500(pucn5s500<=uc5s500c)) length(pindn5s500(pindn5s500<ind5s500c)) length(pccn5s500(pccn5s500<cc5s500c)) length(pind4n5s500(pind4n5s500<ind45s500c)) length(pcc4n5s500(pcc4n5s500<cc45s500c)) length(dq1n5s500(dq1n5s500>dq15s500c)) length(dq4n5s500(dq4n5s500>dq45s500c)) length(tsvqn5s500(tsvqn5s500>tsvq5s500c)) length(BFucn5s500(BFucn5s500>BFuc5s500c)) length(BFind1n5s500(BFind1n5s500>BFind5s500c)) length(BFcc1n5s500(BFcc1n5s500>BFcc5s500c)) length(BFind4n5s500(BFind4n5s500>BFind45s500c)) length(BFcc4n5s500(BFcc4n5s500>BFcc45s500c)) length(BFvqn5s500(BFvqn5s500>BFvq5s500c)) length(BFdq1n5s500(BFdq1n5s500>BFdq15s500c)) length(BFdq14n5s500(BFdq14n5s500>BFdq145s500c)) length(BFdq2n5s500(BFdq2n5s500>BFdq25s500c))]./length(pvqrn5s500)
[length(pucn1s500(pucn1s500<=uc1s500c)) length(pindn1s500(pindn1s500<ind1s500c)) length(pccn1s500(pccn1s500<cc1s500c)) length(pind4n1s500(pind4n1s500<ind41s500c)) length(pcc4n1s500(pcc4n1s500<cc41s500c)) length(dq1n1s500(dq1n1s500>dq11s500c)) length(dq4n1s500(dq4n1s500>dq41s500c)) length(tsvqn1s500(tsvqn1s500>tsvq1s500c)) length(BFucn1s500(BFucn1s500>BFuc1s500c)) length(BFind1n1s500(BFind1n1s500>BFind1s500c)) length(BFcc1n1s500(BFcc1n1s500>BFcc1s500c)) length(BFind4n1s500(BFind4n1s500>BFind41s500c)) length(BFcc4n1s500(BFcc4n1s500>BFcc41s500c)) length(BFvqn1s500(BFvqn1s500>BFvq1s500c)) length(BFdq1n1s500(BFdq1n1s500>BFdq15s500c)) length(BFdq14n1s500(BFdq14n1s500>BFdq141s500c)) length(BFdq2n1s500(BFdq2n1s500>BFdq21s500c))]./length(pvqrn1s500)
%[length(pucn5s500(pucn5s500<=uc5s500c)) length(pindn5s500(pindn5s500<ind5s500c)) length(pccn5s500(pccn5s500<cc5s500c)) length(pind4n5s500(pind4n5s500<ind45s500c)) length(pcc4n5s500(pcc4n5s500<cc45s500c)) length(dq1n5s500(dq1n5s500>dq15s500c)) length(dq4n5s500(dq4n5s500>dq45s500c)) length(tsvqn5s500(tsvqn5s500>tsvq5s500c)) length(BFucn5s500(BFucn5s500>BFuc5s500c)) length(BFind1n5s500(BFind1n5s500>BFind5s500c)) length(BFcc1n5s500(BFcc1n5s500>BFcc5s500c)) length(BFind4n5s500(BFind4n5s500>BFind45s500c)) length(BFcc4n5s500(BFcc4n5s500>BFcc45s500c)) length(BFvqn5s500(BFvqn5s500>BFvq5s500c)) length(BFdq1n5s500(BFdq1n5s500>BFdq15s500c)) length(BFdq14n5s500(BFdq14n5s500>BFdq145s500c))]./length(pvqrn5s500)
%[length(pucn1s500(pucn1s500<=uc1s500c)) length(pindn1s500(pindn1s500<ind1s500c)) length(pccn1s500(pccn1s500<cc1s500c)) length(pind4n1s500(pind4n1s500<ind41s500c)) length(pcc4n1s500(pcc4n1s500<cc41s500c)) length(dq1n1s500(dq1n1s500>dq11s500c)) length(dq4n1s500(dq4n1s500>dq41s500c)) length(tsvqn1s500(tsvqn1s500>tsvq1s500c)) length(BFucn1s500(BFucn1s500>BFuc1s500c)) length(BFind1n1s500(BFind1n1s500>BFind1s500c)) length(BFcc1n1s500(BFcc1n1s500>BFcc1s500c)) length(BFind4n1s500(BFind4n1s500>BFind41s500c)) length(BFcc4n1s500(BFcc4n1s500>BFcc41s500c)) length(BFvqn1s500(BFvqn1s500>BFvq1s500c)) length(BFdq1n1s500(BFdq1n1s500>BFdq11s500c)) length(BFdq14n1s500(BFdq14n1s500>BFdq141s500c))]./length(pvqrn1s500)


disp('size-adjusted power');
[length(pucn5p500(pucn5p500<=uc5s500c)) length(pindn5p500(pindn5p500<ind5s500c)) length(pccn5p500(pccn5p500<cc5s500c)) length(pind4n5p500(pind4n5p500<ind45s500c)) length(pcc4n5p500(pcc4n5p500<cc45s500c)) length(dq1n5p500(dq1n5p500>dq15s500c)) length(dq4n5p500(dq4n5p500>dq45s500c)) length(tsvqn5p500(tsvqn5p500>tsvq5s500c)) length(BFucn5p500(BFucn5p500>BFuc5s500c)) length(BFind1n5p500(BFind1n5p500>BFind5s500c)) length(BFcc1n5p500(BFcc1n5p500>BFcc5s500c)) length(BFind4n5p500(BFind4n5p500>BFind45s500c)) length(BFcc4n5p500(BFcc4n5p500>BFcc45s500c)) length(BFvqn5p500(BFvqn5p500>BFvq5s500c)) length(BFdq1n5p500(BFdq1n5p500>BFdq15s500c)) length(BFdq14n5p500(BFdq14n5p500>BFdq145s500c)) length(BFdq2n5p500(BFdq2n5p500>BFdq25s500c))]./length(pvqrn5p500)
[length(pucn1p500(pucn1p500<=uc1s500c)) length(pindn1p500(pindn1p500<ind1s500c)) length(pccn1p500(pccn1p500<cc1s500c)) length(pind4n1p500(pind4n1p500<ind41s500c)) length(pcc4n1p500(pcc4n1p500<cc41s500c)) length(dq1n1p500(dq1n1p500>dq11s500c)) length(dq4n1p500(dq4n1p500>dq41s500c)) length(tsvqn1p500(tsvqn1p500>tsvq1s500c)) length(BFucn1p500(BFucn1p500>BFuc1s500c)) length(BFind1n1p500(BFind1n1p500>BFind1s500c)) length(BFcc1n1p500(BFcc1n1p500>BFcc1s500c)) length(BFind4n1p500(BFind4n1p500>BFind41s500c)) length(BFcc4n1p500(BFcc4n1p500>BFcc41s500c)) length(BFvqn1p500(BFvqn1p500>BFvq1s500c)) length(BFdq1n1p500(BFdq1n1p500>BFdq11s500c)) length(BFdq14n1p500(BFdq14n1p500>BFdq141s500c)) length(BFdq2n1p500(BFdq2n1p500>BFdq21s500c))]./length(pvqrn1p500)
%[length(pucn5p500(pucn5p500<=uc5s500c)) length(pindn5p500(pindn5p500<ind5s500c)) length(pccn5p500(pccn5p500<cc5s500c)) length(pind4n5p500(pind4n5p500<ind45s500c)) length(pcc4n5p500(pcc4n5p500<cc45s500c)) length(dq1n5p500(dq1n5p500>dq15s500c)) length(dq4n5p500(dq4n5p500>dq45s500c)) length(tsvqn5p500(tsvqn5p500>tsvq5s500c)) length(BFucn5p500(BFucn5p500>BFuc5s500c)) length(BFind1n5p500(BFind1n5p500>BFind5s500c)) length(BFcc1n5p500(BFcc1n5p500>BFcc5s500c)) length(BFind4n5p500(BFind4n5p500>BFind45s500c)) length(BFcc4n5p500(BFcc4n5p500>BFcc45s500c)) length(BFvqn5p500(BFvqn5p500>BFvq5s500c)) length(BFdq1n5p500(BFdq1n5p500>BFdq15s500c)) length(BFdq14n5p500(BFdq14n5p500>BFdq145s500c))]./length(pvqrn5p500)
%[length(pucn1p500(pucn1p500<=uc1s500c)) length(pindn1p500(pindn1p500<ind1s500c)) length(pccn1p500(pccn1p500<cc1s500c)) length(pind4n1p500(pind4n1p500<ind41s500c)) length(pcc4n1p500(pcc4n1p500<cc41s500c)) length(dq1n1p500(dq1n1p500>dq11s500c)) length(dq4n1p500(dq4n1p500>dq41s500c)) length(tsvqn1p500(tsvqn1p500>tsvq1s500c)) length(BFucn1p500(BFucn1p500>BFuc1s500c)) length(BFind1n1p500(BFind1n1p500>BFind1s500c)) length(BFcc1n1p500(BFcc1n1p500>BFcc1s500c)) length(BFind4n1p500(BFind4n1p500>BFind41s500c)) length(BFcc4n1p500(BFcc4n1p500>BFcc41s500c)) length(BFvqn1p500(BFvqn1p500>BFvq1s500c)) length(BFdq1n1p500(BFdq1n1p500>BFdq11s500c)) length(BFdq14n1p500(BFdq14n1p500>BFdq141s500c))]./length(pvqrn1p500)

Ms5500=[hits1s500' pucn5s500' lrucn5s500' pindn5s500' lrindn5s500' pccn5s500' lrccn5s500' pind4n5s500' lrind4n5s500' pcc4n5s500' lrcc4n5s500' dq1n5s500' pdq1n5s500' dq4n5s500' pdq4n5s500' tsvqn5s500' pvqrn5s500' BFucn5s500' BFind1n5s500' BFcc1n5s500' BFind4n5s500' BFcc4n5s500' BFvqn5s500' BFdq1n5s500' BFdq14n5s500' hits5s500' BFdq2n5s500']; 
%Ms5500=[hits5s500' pucn5s500' lrucn5s500' pindn5s500' lrindn5s500' pccn5s500' lrccn5s500' pind4n5s500' lrind4n5s500' pcc4n5s500' lrcc4n5s500' dq1n5s500' pdq1n5s500' dq4n5s500' pdq4n5s500' tsvqn5s500' pvqrn5s500' BFucn5s500' BFind1n5s500' BFcc1n5s500' BFind4n5s500' BFcc4n5s500' BFvqn5s500']; 
csvwrite('HS500s_5final2.csv',Ms5500)
Ms1500=[hits1s500' pucn1s500' lrucn1s500' pindn1s500' lrindn1s500' pccn1s500' lrccn1s500' pind4n1s500' lrind4n1s500' pcc4n1s500' lrcc4n1s500' dq1n1s500' pdq1n1s500' dq4n1s500' pdq4n1s500' tsvqn1s500' pvqrn1s500' BFucn1s500' BFind1n1s500' BFcc1n1s500' BFind4n1s500' BFcc4n1s500' BFvqn1s500' BFdq1n1s500' BFdq14n1s500' hits1s500' BFdq2n1s500']; 
%Ms1500=[hits1s500' pucn1s500' lrucn1s500' pindn1s500' lrindn1s500' pccn1s500' lrccn1s500' pind4n1s500' lrind4n1s500' pcc4n1s500' lrcc4n1s500' dq1n1s500' pdq1n1s500' dq4n1s500' pdq4n1s500' tsvqn1s500' pvqrn1s500' BFucn1s500' BFind1n1s500' BFcc1n1s500' BFind4n1s500' BFcc4n1s500' BFvqn1s500']; 
csvwrite('HS500s_1final2.csv',Ms1500)
Mp5500=[hits5p500' pucn5p500' lrucn5p500' pindn5p500' lrindn5p500' pccn5p500' lrccn5p500' pind4n5p500' lrind4n5p500' pcc4n5p500' lrcc4n5p500' dq1n5p500' pdq1n5p500' dq4n5p500' pdq4n5p500' tsvqn5p500' pvqrn5p500' BFucn5p500' BFind1n5p500' BFcc1n5p500' BFind4n5p500' BFcc4n5p500' BFvqn5p500' BFdq1n5p500' BFdq14n5p500' hits5p500' BFdq2n5p500']; 
%Mp5500=[hits5p500' pucn5p500' lrucn5p500' pindn5p500' lrindn5p500' pccn5p500' lrccn5p500' pind4n5p500' lrind4n5p500' pcc4n5p500' lrcc4n5p500' dq1n5p500' pdq1n5p500' dq4n5p500' pdq4n5p500' tsvqn5p500' pvqrn5p500' BFucn5p500' BFind1n5p500' BFcc1n5p500' BFind4n5p500' BFcc4n5p500' BFvqn5p500']; 
csvwrite('HS500p_5final2.csv',Mp5500)
Mp1500=[hits1p500' pucn1p500' lrucn1p500' pindn1p500' lrindn1p500' pccn1p500' lrccn1p500' pind4n1p500' lrind4n1p500' pcc4n1p500' lrcc4n1p500' dq1n1p500' pdq1n1p500' dq4n1p500' pdq4n1p500' tsvqn1p500' pvqrn1p500' BFucn1p500' BFind1n1p500' BFcc1n1p500' BFind4n1p500' BFcc4n1p500' BFvqn1p500' BFdq1n1p500' BFdq14n1p500' BFdq2n1p500']; 
%Mp1500=[hits1p500' pucn1p500' lrucn1p500' pindn1p500' lrindn1p500' pccn1p500' lrccn1p500' pind4n1p500' lrind4n1p500' pcc4n1p500' lrcc4n1p500' dq1n1p500' pdq1n1p500' dq4n1p500' pdq4n1p500' tsvqn1p500' pvqrn1p500' BFucn1p500' BFind1n1p500' BFcc1n1p500' BFind4n1p500' BFcc4n1p500' BFvqn1p500']; 
csvwrite('HS500p_1final2.csv',Mp1500)

n=250;
for j=1:20000
j
[y, sig, VaR5, VaR1, VaR5HS, VaR1HS]=gengarchnHS2(n,a0,a1,b1);
  hits5s250(j)=sum(y<VaR5);hits5p250(j)=sum(y<VaR5HS);
  hits1s250(j)=sum(y<VaR1);hits1p250(j)=sum(y<VaR1HS);
  [pucn1s250(j),lrucn1s250(j),pccn1s250(j),lrccn1s250(j),pindn1s250(j),lrindn1s250(j),pcc4n1s250(j),lrcc4n1s250(j),pind4n1s250(j),lrind4n1s250(j),pdq1n1s250(j),dq1n1s250(j),pdq4n1s250(j),dq4n1s250(j),tsvqn1s250(j),pvqrn1s250(j)]=runucccdqvqtests(y,VaR1,0.01);
  [pucn5s250(j),lrucn5s250(j),pccn5s250(j),lrccn5s250(j),pindn5s250(j),lrindn5s250(j),pcc4n5s250(j),lrcc4n5s250(j),pind4n5s250(j),lrind4n5s250(j),pdq1n5s250(j),dq1n5s250(j),pdq4n5s250(j),dq4n5s250(j),tsvqn5s250(j),pvqrn5s250(j)]=runucccdqvqtests(y,VaR5,0.05);
  [pucn1p250(j),lrucn1p250(j),pccn1p250(j),lrccn1p250(j),pindn1p250(j),lrindn1p250(j),pcc4n1p250(j),lrcc4n1p250(j),pind4n1p250(j),lrind4n1p250(j),pdq1n1p250(j),dq1n1p250(j),pdq4n1p250(j),dq4n1p250(j),tsvqn1p250(j),pvqrn1p250(j)]=runucccdqvqtests(y,VaR1HS,0.01);
  [pucn5p250(j),lrucn5p250(j),pccn5p250(j),lrccn5p250(j),pindn5p250(j),lrindn5p250(j),pcc4n5p250(j),lrcc4n5p250(j),pind4n5p250(j),lrind4n5p250(j),pdq1n5p250(j),dq1n5p250(j),pdq4n5p250(j),dq4n5p250(j),tsvqn5p250(j),pvqrn5p250(j)]=runucccdqvqtests(y,VaR5HS,0.05);
  [BFucn1s250(j),BFcc1n1s250(j),BFind1n1s250(j),BFcc4n1s250(j),BFind4n1s250(j),BFvqn1s250(j)]=runBFtests(y,VaR1,0.01);
  [BFucn5s250(j),BFcc1n5s250(j),BFind1n5s250(j),BFcc4n5s250(j),BFind4n5s250(j),BFvqn5s250(j)]=runBFtests(y,VaR5,0.05);
 [BFucn1p250(j),BFcc1n1p250(j),BFind1n1p250(j),BFcc4n1p250(j),BFind4n1p250(j),BFvqn1p250(j)]=runBFtests(y,VaR1HS,0.01);
 [BFucn5p250(j),BFcc1n5p250(j),BFind1n5p250(j),BFcc4n5p250(j),BFind4n5p250(j),BFvqn5p250(j)]=runBFtests(y,VaR5HS,0.05);
  BFdq1n1s250(j) = 1/bfdqmlik3(y,VaR1,0.01);BFdq1n5s250(j)=1/bfdqmlik3(y,VaR5,0.05);
  BFdq1n1p250(j) = 1/bfdqmlik3(y,VaR1HS,0.01);BFdq1n5p250(j)=1/bfdqmlik3(y,VaR5HS,0.05);
  BFdq14n1s250(j) = 1/bfdqmlik4(y,VaR1,0.01);BFdq14n5s250(j)=1/bfdqmlik4(y,VaR5,0.05);
  BFdq14n1p250(j) = 1/bfdqmlik4(y,VaR1HS,0.01);BFdq14n5p250(j)=1/bfdqmlik4(y,VaR5HS,0.05);
   BFdq2n1s250(j) = 1/bfdq2mlik2(y,VaR1,0.01);BFdq2n5s250(j)=1/bfdq2mlik2(y,VaR5,0.05);
   BFdq2n1p250(j) = 1/bfdq2mlik2(y,VaR1HS,0.01);BFdq2n5p250(j)=1/bfdq2mlik2(y,VaR5HS,0.05);
  
  if rem(j,100)==0
     j
     [pucn1s250(j),lrucn1s250(j),pccn1s250(j),lrccn1s250(j),pindn1s250(j),lrindn1s250(j),pcc4n1s250(j),lrcc4n1s250(j),pind4n1s250(j),lrind4n1s250(j),pdq1n1s250(j),dq1n1s250(j),pdq4n1s250(j),dq4n1s250(j),tsvqn1s250(j),pvqrn1s250(j);
      pucn1p250(j),lrucn1p250(j),pccn1p250(j),lrccn1p250(j),pindn1p250(j),lrindn1p250(j),pcc4n1p250(j),lrcc4n1p250(j),pind4n1p250(j),lrind4n1p250(j),pdq1n1p250(j),dq1n1p250(j),pdq4n1p250(j),dq4n1p250(j),tsvqn1p250(j),pvqrn1p250(j)]
     %[2*log(BFdq2n1s250(j)) 2*log(BFdq2n5s250(j)) 2*log(BFdq14n1s250(j)) 2*log(BFdq14n5s250(j));
     %  2*log(BFdq2n1p250(j)) 2*log(BFdq2n5p250(j)) 2*log(BFdq14n1p250(j)) 2*log(BFdq14n5p250(j))]  
  end   
  
end  

'size'
[length(pucn5s250(pucn5s250<0.05)) length(pindn5s250(pindn5s250<0.05)) length(pccn5s250(pccn5s250<0.05)) length(pind4n5s250(pind4n5s250<0.05)) length(pcc4n5s250(pcc4n5s250<0.05)) length(pdq1n5s250(pdq1n5s250<0.05)) length(pdq4n5s250(pdq4n5s250<0.05)) length(pvqrn5s250(pvqrn5s250<0.05)) length(BFucn5s250(BFucn5s250>1)) length(BFind1n5s250(BFind1n5s250>1)) length(BFcc1n5s250(BFcc1n5s250>1)) length(BFind4n5s250(BFind4n5s250>1)) length(BFcc4n5s250(BFcc4n5s250>1)) length(BFvqn5s250(BFvqn5s250>1)) length(BFdq1n5s250(BFdq1n5s250>1)) length(BFdq14n5s250(BFdq14n5s250>1)) length(BFdq2n5s250(BFdq2n5s250>1))]./length(pvqrn5s250)
[length(pucn1s250(pucn1s250<0.05)) length(pindn1s250(pindn1s250<0.05)) length(pccn1s250(pccn1s250<0.05)) length(pind4n1s250(pind4n1s250<0.05)) length(pcc4n1s250(pcc4n1s250<0.05)) length(pdq1n1s250(pdq1n1s250<0.05)) length(pdq4n1s250(pdq4n1s250<0.05)) length(pvqrn1s250(pvqrn1s250<0.05)) length(BFucn1s250(BFucn1s250>1)) length(BFind1n1s250(BFind1n1s250>1)) length(BFcc1n1s250(BFcc1n1s250>1))  length(BFind4n1s250(BFind4n1s250>1)) length(BFcc4n1s250(BFcc4n1s250>1)) length(BFvqn1s250(BFvqn1s250>1)) length(BFdq1n1s250(BFdq1n1s250>1)) length(BFdq14n1s250(BFdq14n1s250>1)) length(BFdq2n1s250(BFdq2n1s250>1))]./length(pvqrn1s250)
%[length(pucn5s250(pucn5s250<0.05)) length(pindn5s250(pindn5s250<0.05)) length(pccn5s250(pccn5s250<0.05)) length(pind4n5s250(pind4n5s250<0.05)) length(pcc4n5s250(pcc4n5s250<0.05)) length(pdq1n5s250(pdq1n5s250<0.05)) length(pdq4n5s250(pdq4n5s250<0.05)) length(pvqrn5s250(pvqrn5s250<0.05)) length(BFucn5s250(BFucn5s250>1)) length(BFind1n5s250(BFind1n5s250>1)) length(BFcc1n5s250(BFcc1n5s250>1)) length(BFind4n5s250(BFind4n5s250>1)) length(BFcc4n5s250(BFcc4n5s250>1)) length(BFvqn5s250(BFvqn5s250>1)) length(BFdq1n5s250(BFdq1n5s250>1)) length(BFdq14n5s250(BFdq14n5s250>1))]./length(pvqrn5s250)
%[length(pucn1s250(pucn1s250<0.05)) length(pindn1s250(pindn1s250<0.05)) length(pccn1s250(pccn1s250<0.05)) length(pind4n1s250(pind4n1s250<0.05)) length(pcc4n1s250(pcc4n1s250<0.05)) length(pdq1n1s250(pdq1n1s250<0.05)) length(pdq4n1s250(pdq4n1s250<0.05)) length(pvqrn1s250(pvqrn1s250<0.05)) length(BFucn1s250(BFucn1s250>1)) length(BFind1n1s250(BFind1n1s250>1)) length(BFcc1n1s250(BFcc1n1s250>1))  length(BFind4n1s250(BFind4n1s250>1)) length(BFcc4n1s250(BFcc4n1s250>1)) length(BFvqn1s250(BFvqn1s250>1)) length(BFdq1n1s250(BFdq1n1s250>1)) length(BFdq14n1s250(BFdq14n1s250>1))]./length(pvqrn1s250)


'power'
[length(pucn5p250(pucn5p250<0.05)) length(pindn5p250(pindn5p250<0.05)) length(pccn5p250(pccn5p250<0.05)) length(pind4n5p250(pind4n5p250<0.05)) length(pcc4n5p250(pcc4n5p250<0.05)) length(pdq1n5p250(pdq1n5p250<0.05)) length(pdq4n5p250(pdq4n5p250<0.05)) length(pvqrn5p250(pvqrn5p250<0.05)) length(BFucn5p250(BFucn5p250>1)) length(BFind1n5p250(BFind1n5p250>1)) length(BFcc1n5p250(BFcc1n5p250>1)) length(BFind4n5p250(BFind4n5p250>1)) length(BFcc4n5p250(BFcc4n5p250>1)) length(BFvqn5p250(BFvqn5p250>1)) length(BFdq1n5p250(BFdq1n5p250>1)) length(BFdq14n5p250(BFdq14n5p250>1)) length(BFdq2n5p250(BFdq2n5p250>1))]./length(pvqrn5p250)
[length(pucn1p250(pucn1p250<0.05)) length(pindn1p250(pindn1p250<0.05)) length(pccn1p250(pccn1p250<0.05)) length(pind4n1p250(pind4n1p250<0.05)) length(pcc4n1p250(pcc4n1p250<0.05)) length(pdq1n1p250(pdq1n1p250<0.05)) length(pdq4n1p250(pdq4n1p250<0.05)) length(pvqrn1p250(pvqrn1p250<0.05)) length(BFucn1p250(BFucn1p250>1)) length(BFind1n1p250(BFind1n1p250>1)) length(BFcc1n1p250(BFcc1n1p250>1)) length(BFind4n1p250(BFind4n1p250>1)) length(BFcc4n1p250(BFcc4n1p250>1)) length(BFvqn1p250(BFvqn1p250>1)) length(BFdq1n1p250(BFdq1n1p250>1)) length(BFdq14n1p250(BFdq14n5p250>1)) length(BFdq2n1p250(BFdq2n1p250>1))]./length(pvqrn1p250)
%[length(pucn5p250(pucn5p250<0.05)) length(pindn5p250(pindn5p250<0.05)) length(pccn5p250(pccn5p250<0.05)) length(pind4n5p250(pind4n5p250<0.05)) length(pcc4n5p250(pcc4n5p250<0.05)) length(pdq1n5p250(pdq1n5p250<0.05)) length(pdq4n5p250(pdq4n5p250<0.05)) length(pvqrn5p250(pvqrn5p250<0.05)) length(BFucn5p250(BFucn5p250>1)) length(BFind1n5p250(BFind1n5p250>1)) length(BFcc1n5p250(BFcc1n5p250>1)) length(BFind4n5p250(BFind4n5p250>1)) length(BFcc4n5p250(BFcc4n5p250>1)) length(BFvqn5p250(BFvqn5p250>1)) length(BFdq1n5p250(BFdq1n5p250>1)) length(BFdq14n5p250(BFdq14n5p250>1))]./length(pvqrn5p250)
%[length(pucn1p250(pucn1p250<0.05)) length(pindn1p250(pindn1p250<0.05)) length(pccn1p250(pccn1p250<0.05)) length(pind4n1p250(pind4n1p250<0.05)) length(pcc4n1p250(pcc4n1p250<0.05)) length(pdq1n1p250(pdq1n1p250<0.05)) length(pdq4n1p250(pdq4n1p250<0.05)) length(pvqrn1p250(pvqrn1p250<0.05)) length(BFucn1p250(BFucn1p250>1)) length(BFind1n1p250(BFind1n1p250>1)) length(BFcc1n1p250(BFcc1n1p250>1)) length(BFind4n1p250(BFind4n1p250>1)) length(BFcc4n1p250(BFcc4n1p250>1)) length(BFvqn1p250(BFvqn1p250>1)) length(BFdq1n1p250(BFdq1n1p250>1)) length(BFdq14n1p250(BFdq14n1p250>1))]./length(pvqrn1p250)

BFuc5s250c=prctile(BFucn5s250,95);BFind5s250c=prctile(BFind1n5s250,95);
BFcc5s250c=prctile(BFcc1n5s250,95);BFcc45s250c=prctile(BFcc4n5s250,95);
BFind45s250c=prctile(BFind4n5s250,95);BFvq5s250c=prctile(BFvqn5s250,95);
BFdq15s250c=prctile(BFdq1n5s250,95);BFdq145s250c=prctile(BFdq14n5s250,95);BFdq25s250c=prctile(BFdq2n5s250,95);
uc5s250c=prctile(pucn5s250,5);ind5s250c=prctile(pindn5s250,5);ind45s250c=prctile(pind4n5s250,5);
cc5s250c=prctile(pccn5s250,5);cc45s250c=prctile(pcc4n5s250,5);
dq15s250c=prctile(dq1n5s250,95);dq45s250c=prctile(dq4n5s250,95);tsvq5s250c=prctile(tsvqn5s250,95);

BFuc1s250c=prctile(BFucn1s250,95);BFind1s250c=prctile(BFind1n1s250,95);
BFcc1s250c=prctile(BFcc1n1s250,95);BFcc41s250c=prctile(BFcc4n1s250,95);
BFind41s250c=prctile(BFind4n1s250,95);BFvq1s250c=prctile(BFvqn1s250,95);
BFdq11s250c=prctile(BFdq1n1s250,95);BFdq141s250c=prctile(BFdq14n1s250,95);BFdq21s250c=prctile(BFdq2n1s250,95);
uc1s250c=prctile(pucn1s250,5);ind1s250c=prctile(pindn1s250,5);ind41s250c=prctile(pind4n1s250,5);
cc1s250c=prctile(pccn1s250,5);cc41s250c=prctile(pcc4n1s250,5);
dq11s250c=prctile(dq1n1s250,95);dq41s250c=prctile(dq4n1s250,95);tsvq1s250c=prctile(tsvqn1s250,95);

'corrected size'
[length(pucn5s250(pucn5s250<=uc5s250c)) length(pindn5s250(pindn5s250<ind5s250c)) length(pccn5s250(pccn5s250<cc5s250c)) length(pind4n5s250(pind4n5s250<ind45s250c)) length(pcc4n5s250(pcc4n5s250<cc45s250c)) length(dq1n5s250(dq1n5s250>dq15s250c)) length(dq4n5s250(dq4n5s250>dq45s250c)) length(tsvqn5s250(tsvqn5s250>tsvq5s250c)) length(BFucn5s250(BFucn5s250>BFuc5s250c)) length(BFind1n5s250(BFind1n5s250>BFind5s250c)) length(BFcc1n5s250(BFcc1n5s250>BFcc5s250c)) length(BFind4n5s250(BFind4n5s250>BFind45s250c)) length(BFcc4n5s250(BFcc4n5s250>BFcc45s250c)) length(BFvqn5s250(BFvqn5s250>BFvq5s250c)) length(BFdq1n5s250(BFdq1n5s250>BFdq15s250c)) length(BFdq14n5s250(BFdq14n5s250>BFdq145s250c)) length(BFdq2n5s250(BFdq2n5s250>BFdq25s250c))]./length(pvqrn5s250)
[length(pucn1s250(pucn1s250<=uc1s250c)) length(pindn1s250(pindn1s250<ind1s250c)) length(pccn1s250(pccn1s250<cc1s250c)) length(pind4n1s250(pind4n1s250<ind41s250c)) length(pcc4n1s250(pcc4n1s250<cc41s250c)) length(dq1n1s250(dq1n1s250>dq11s250c)) length(dq4n1s250(dq4n1s250>dq41s250c)) length(tsvqn1s250(tsvqn1s250>tsvq1s250c)) length(BFucn1s250(BFucn1s250>BFuc1s250c)) length(BFind1n1s250(BFind1n1s250>BFind1s250c)) length(BFcc1n1s250(BFcc1n1s250>BFcc1s250c)) length(BFind4n1s250(BFind4n1s250>BFind41s250c)) length(BFcc4n1s250(BFcc4n1s250>BFcc41s250c)) length(BFvqn1s250(BFvqn1s250>BFvq1s250c)) length(BFdq1n1s250(BFdq1n1s250>BFdq11s250c)) length(BFdq14n1s250(BFdq14n1s250>BFdq141s250c)) length(BFdq2n1s250(BFdq2n1s250>BFdq21s250c))]./length(pvqrn1s250)
%[length(pucn5s250(pucn5s250<=uc5s250c)) length(pindn5s250(pindn5s250<ind5s250c)) length(pccn5s250(pccn5s250<cc5s250c)) length(pind4n5s250(pind4n5s250<ind45s250c)) length(pcc4n5s250(pcc4n5s250<cc45s250c)) length(dq1n5s250(dq1n5s250>dq15s250c)) length(dq4n5s250(dq4n5s250>dq45s250c)) length(tsvqn5s250(tsvqn5s250>tsvq5s250c)) length(BFucn5s250(BFucn5s250>BFuc5s250c)) length(BFind1n5s250(BFind1n5s250>BFind5s250c)) length(BFcc1n5s250(BFcc1n5s250>BFcc5s250c)) length(BFind4n5s250(BFind4n5s250>BFind45s250c)) length(BFcc4n5s250(BFcc4n5s250>BFcc45s250c)) length(BFvqn5s250(BFvqn5s250>BFvq5s250c)) length(BFdq1n5s250(BFdq1n5s250>BFdq15s250c)) length(BFdq14n5s250(BFdq14n5s250>BFdq145s250c))]./length(pvqrn5s250)
%[length(pucn1s250(pucn1s250<=uc1s250c)) length(pindn1s250(pindn1s250<ind1s250c)) length(pccn1s250(pccn1s250<cc1s250c)) length(pind4n1s250(pind4n1s250<ind41s250c)) length(pcc4n1s250(pcc4n1s250<cc41s250c)) length(dq1n1s250(dq1n1s250>dq11s250c)) length(dq4n1s250(dq4n1s250>dq41s250c)) length(tsvqn1s250(tsvqn1s250>tsvq1s250c)) length(BFucn1s250(BFucn1s250>BFuc1s250c)) length(BFind1n1s250(BFind1n1s250>BFind1s250c)) length(BFcc1n1s250(BFcc1n1s250>BFcc1s250c)) length(BFind4n1s250(BFind4n1s250>BFind41s250c)) length(BFcc4n1s250(BFcc4n1s250>BFcc41s250c)) length(BFvqn1s250(BFvqn1s250>BFvq1s250c)) length(BFdq1n1s250(BFdq1n1s250>BFdq11s250c)) length(BFdq14n1s250(BFdq14n1s250>BFdq141s250c))]./length(pvqrn1s250)


'size-adjusted power'
[length(pucn5p250(pucn5p250<=uc5s250c)) length(pindn5p250(pindn5p250<ind5s250c)) length(pccn5p250(pccn5p250<cc5s250c)) length(pind4n5p250(pind4n5p250<ind45s250c)) length(pcc4n5p250(pcc4n5p250<cc45s250c)) length(dq1n5p250(dq1n5p250>dq15s250c)) length(dq4n5p250(dq4n5p250>dq45s250c)) length(tsvqn5p250(tsvqn5p250>tsvq5s250c)) length(BFucn5p250(BFucn5p250>BFuc5s250c)) length(BFind1n5p250(BFind1n5p250>BFind5s250c)) length(BFcc1n5p250(BFcc1n5p250>BFcc5s250c)) length(BFind4n5p250(BFind4n5p250>BFind45s250c)) length(BFcc4n5p250(BFcc4n5p250>BFcc45s250c)) length(BFvqn5p250(BFvqn5p250>BFvq5s250c)) length(BFdq1n5p250(BFdq1n5p250>BFdq15s250c)) length(BFdq14n5p250(BFdq14n5p250>BFdq145s250c)) length(BFdq2n5p250(BFdq2n5p250>BFdq25s250c))]./length(pvqrn5p250)
[length(pucn1p250(pucn1p250<=uc1s250c)) length(pindn1p250(pindn1p250<ind1s250c)) length(pccn1p250(pccn1p250<cc1s250c)) length(pind4n1p250(pind4n1p250<ind41s250c)) length(pcc4n1p250(pcc4n1p250<cc41s250c)) length(dq1n1p250(dq1n1p250>dq11s250c)) length(dq4n1p250(dq4n1p250>dq41s250c)) length(tsvqn1p250(tsvqn1p250>tsvq1s250c)) length(BFucn1p250(BFucn1p250>BFuc1s250c)) length(BFind1n1p250(BFind1n1p250>BFind1s250c)) length(BFcc1n1p250(BFcc1n1p250>BFcc1s250c)) length(BFind4n1p250(BFind4n1p250>BFind41s250c)) length(BFcc4n1p250(BFcc4n1p250>BFcc41s250c)) length(BFvqn1p250(BFvqn1p250>BFvq1s250c)) length(BFdq1n1p250(BFdq1n1p250>BFdq11s250c)) length(BFdq14n1p250(BFdq14n1p250>BFdq141s250c)) length(BFdq2n1p250(BFdq2n1p250>BFdq21s250c))]./length(pvqrn1p250)
%[length(pucn5p250(pucn5p250<=uc5s250c)) length(pindn5p250(pindn5p250<ind5s250c)) length(pccn5p250(pccn5p250<cc5s250c)) length(pind4n5p250(pind4n5p250<ind45s250c)) length(pcc4n5p250(pcc4n5p250<cc45s250c)) length(dq1n5p250(dq1n5p250>dq15s250c)) length(dq4n5p250(dq4n5p250>dq45s250c)) length(tsvqn5p250(tsvqn5p250>tsvq5s250c)) length(BFucn5p250(BFucn5p250>BFuc5s250c)) length(BFind1n5p250(BFind1n5p250>BFind5s250c)) length(BFcc1n5p250(BFcc1n5p250>BFcc5s250c)) length(BFind4n5p250(BFind4n5p250>BFind45s250c)) length(BFcc4n5p250(BFcc4n5p250>BFcc45s250c)) length(BFvqn5p250(BFvqn5p250>BFvq5s250c)) length(BFdq1n5p250(BFdq1n5p250>BFdq15s250c)) length(BFdq14n5p250(BFdq14n5p250>BFdq145s250c))]./length(pvqrn5p250)
%[length(pucn1p250(pucn1p250<=uc1s250c)) length(pindn1p250(pindn1p250<ind1s250c)) length(pccn1p250(pccn1p250<cc1s250c)) length(pind4n1p250(pind4n1p250<ind41s250c)) length(pcc4n1p250(pcc4n1p250<cc41s250c)) length(dq1n1p250(dq1n1p250>dq11s250c)) length(dq4n1p250(dq4n1p250>dq41s250c)) length(tsvqn1p250(tsvqn1p250>tsvq1s250c)) length(BFucn1p250(BFucn1p250>BFuc1s250c)) length(BFind1n1p250(BFind1n1p250>BFind1s250c)) length(BFcc1n1p250(BFcc1n1p250>BFcc1s250c)) length(BFind4n1p250(BFind4n1p250>BFind41s250c)) length(BFcc4n1p250(BFcc4n1p250>BFcc41s250c)) length(BFvqn1p250(BFvqn1p250>BFvq1s250c)) length(BFdq1n1p250(BFdq1n1p250>BFdq11s250c)) length(BFdq14n1p250(BFdq14n1p250>BFdq141s250c))]./length(pvqrn1p250)

Ms5250=[hits5s250' pucn5s250' lrucn5s250' pindn5s250' lrindn5s250' pccn5s250' lrccn5s250' pind4n5s250' lrind4n5s250' pcc4n5s250' lrcc4n5s250' dq1n5s250' pdq1n5s250' dq4n5s250' pdq4n5s250' tsvqn5s250' pvqrn5s250' BFucn5s250' BFind1n5s250' BFcc1n5s250' BFind4n5s250' BFcc4n5s250' BFvqn5s250' BFdq1n5s250' BFdq14n5s250' BFdq2n5s250']; 
%Ms5250=[pucn5s250' lrucn5s250' pindn5s250' lrindn5s250' pccn5s250' lrccn5s250' pind4n5s250' lrind4n5s250' pcc4n5s250' lrcc4n5s250' dq1n5s250' pdq1n5s250' dq4n5s250' pdq4n5s250' tsvqn5s250' pvqrn5s250' BFucn5s250' BFind1n5s250' BFcc1n5s250' BFind4n5s250' BFcc4n5s250' BFvqn5s250']; 
csvwrite('HS250s_5final4.csv',Ms5250)
Ms1250=[hits1s250' pucn1s250' lrucn1s250' pindn1s250' lrindn1s250' pccn1s250' lrccn1s250' pind4n1s250' lrind4n1s250' pcc4n1s250' lrcc4n1s250' dq1n1s250' pdq1n1s250' dq4n1s250' pdq4n1s250' tsvqn1s250' pvqrn1s250' BFucn1s250' BFind1n1s250' BFcc1n1s250' BFind4n1s250' BFcc4n1s250' BFvqn1s250' BFdq1n1s250' BFdq14n1s250' BFdq2n1s250']; 
%Ms1250=[hits1s250' pucn1s250' lrucn1s250' pindn1s250' lrindn1s250' pccn1s250' lrccn1s250' pind4n1s250' lrind4n1s250' pcc4n1s250' lrcc4n1s250' dq1n1s250' pdq1n1s250' dq4n1s250' pdq4n1s250' tsvqn1s250' pvqrn1s250' BFucn1s250' BFind1n1s250' BFcc1n1s250' BFind4n1s250' BFcc4n1s250' BFvqn1s250']; 
csvwrite('HS250s_1final4.csv',Ms1250)
Mp5250=[hits5p250' pucn5p250' lrucn5p250' pindn5p250' lrindn5p250' pccn5p250' lrccn5p250' pind4n5p250' lrind4n5p250' pcc4n5p250' lrcc4n5p250' dq1n5p250' pdq1n5p250' dq4n5p250' pdq4n5p250' tsvqn5p250' pvqrn5p250' BFucn5p250' BFind1n5p250' BFcc1n5p250' BFind4n5p250' BFcc4n5p250' BFvqn5p250' BFdq1n5p250' BFdq14n5p250' BFdq2n5p250']; 
%Mp5250=[hits5p250' pucn5p250' lrucn5p250' pindn5p250' lrindn5p250' pccn5p250' lrccn5p250' pind4n5p250' lrind4n5p250' pcc4n5p250' lrcc4n5p250' dq1n5p250' pdq1n5p250' dq4n5p250' pdq4n5p250' tsvqn5p250' pvqrn5p250' BFucn5p250' BFind1n5p250' BFcc1n5p250' BFind4n5p250' BFcc4n5p250' BFvqn5p250']; 
csvwrite('HS250p_5final4.csv',Mp5250)
Mp1250=[hits1p250' pucn1p250' lrucn1p250' pindn1p250' lrindn1p250' pccn1p250' lrccn1p250' pind4n1p250' lrind4n1p250' pcc4n1p250' lrcc4n1p250' dq1n1p250' pdq1n1p250' dq4n1p250' pdq4n1p250' tsvqn1p250' pvqrn1p250' BFucn1p250' BFind1n1p250' BFcc1n1p250' BFind4n1p250' BFcc4n1p250' BFvqn1p250' BFdq1n1p250' BFdq14n1p250' BFdq2n1p250']; 
%Mp1250=[hits1p250' pucn1p250' lrucn1p250' pindn1p250' lrindn1p250' pccn1p250' lrccn1p250' pind4n1p250' lrind4n1p250' pcc4n1p250' lrcc4n1p250' dq1n1p250' pdq1n1p250' dq4n1p250' pdq4n1p250' tsvqn1p250' pvqrn1p250' BFucn1p250' BFind1n1p250' BFcc1n1p250' BFind4n1p250' BFcc4n1p250' BFvqn1p250']; 
csvwrite('HS250p_1final4.csv',Mp1250)

%empirical 5% critical values or 5% sig. levels

[uc5s250c; uc5s250c; uc5s250c; uc5s250c]
[ind5s250c; ind5s250c; ind5s250c; ind5s250c]
[cc5s250c; cc5s250c; cc5s250c; cc5s250c]
[ind45s250c; ind45s250c; ind45s250c; ind45s250c]
[cc45s250c; cc45s250c; cc45s250c; cc45s250c]
[dq15s250c; dq15s250c; dq15s250c; dq15s250c]
[dq45s250c; dq45s250c; dq45s250c; dq45s250c]
[tsvq5s250c; tsvq5s250c; tsvq5s250c; tsvq5s250c]
[BFuc5s250c; BFuc5s250c; BFuc5s250c; BFuc5s250c]
[BFind5s250c; BFind5s250c; BFind5s250c; BFind5s250c]
[BFcc5s250c; BFcc5s250c; BFcc5s250c; BFcc5s250c]
[BFind45s250c; BFind45s250c; BFind45s250c; BFind45s250c]
[BFcc45s250c; BFcc45s250c; BFcc45s250c; BFcc45s250c]
[BFvq5s250c; BFvq5s250c; BFvq5s250c; BFvq5s250c]

[uc1s250c; uc1s250c; uc1s250c; uc1s250c]
[ind1s250c; ind1s250c; ind1s250c; ind1s250c]
[cc1s250c; cc1s250c; cc1s250c; cc1s250c]
[ind41s250c; ind41s250c; ind41s250c; ind41s250c]
[cc41s250c; cc41s250c; cc41s250c; cc41s250c]
[dq11s250c; dq11s250c; dq11s250c; dq11s250c]
[dq41s250c; dq41s250c; dq41s250c; dq41s250c]
[tsvq1s250c; tsvq1s250c; tsvq1s250c; tsvq1s250c]
[BFuc1s250c; BFuc1s250c; BFuc1s250c; BFuc1s250c]
[BFind1s250c; BFind1s250c; BFind1s250c; BFind1s250c]
[BFcc1s250c; BFcc1s250c; BFcc1s250c; BFcc1s250c]
[BFind41s250c; BFind41s250c; BFind41s250c; BFind41s250c]
[BFcc41s250c; BFcc41s250c; BFcc41s250c; BFcc41s250c]
[BFvq1s250c; BFvq1s250c; BFvq1s250c; BFvq1s250c]
