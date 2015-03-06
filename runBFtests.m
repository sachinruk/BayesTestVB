function [BFuc,BFcc1,BFind1,BFcc4,BFind4,BFvq]=runBFtests(y,f,p)

hits=(y<f);
%you sure you want to ignore the BFuc step from here?
[~,BFcc1,BFind1]=ucccBFtests(hits,p,1); 
[BFuc,BFcc4,BFind4]=ucccBFtests(hits,p,4);
BFvq = 1/bfmlik2_s(y,f,p);