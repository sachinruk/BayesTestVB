function  bf=bfmlik2_s(y,fy,alp)
persistent x0 x1;

    
n=length(y);

% g0=-0.99:0.01:0.99;ng=length(g0);
% g1=g0;
% for i=1:ng
%     for j=1:ng
%        lmlik(i,j)=qrmlik(g0(i),g1(j),y,fy,alp,n);
%     end
% end
%if sum(y<fy)>0  
%mlik=exp(lmlik-max(max(lmlik)));
% A=[eye(2);-1*eye(2)];B=[0.99 0.99 0.99 0.99]';
% X0=[-0.1 0.5]';
% [~,xm]=fmincon(@(g2) qrmlik2(g2,y,fy,alp,n),X0,A,B);
%marginal likelihood, double integral

% xm=-xm;
dx=0.01; dy=dx;
if isempty(x0)
    lim=0.999;
    [X0, X1]=meshgrid(-lim:dx:lim,-lim:dy:lim);
    x0=X0(:); x1=X1(:); 
end
%     n2=length(x0);
a=qrmlik1_s(x0',x1',y,fy,alp,n); 
lmlik2 =logsumexp(a')+log(dx)+log(dy);
% lmlik2=integral2(@(x0,x1)qrmlik1_s(x0,x1,y,fy,alp,n),-0.999,0.999,-0.999,0.999);
% figure; surf(X0,X1,reshape(a,size(X0)))
% mlik2 = dblquad(@(x0,x1) qrmlik1(x0,x1,y,fy,alp,n,xm),-0.99,0.99,-0.99,0.99);

disp(lmlik2);
%Bayes factor

%likelihood  b0=0, b1=1;

u=y-fy;
llik= -n*log(sum(u.*(alp-(u<0))))+gammaln(n)+n*log(alp*(1-alp));
llik1=llik;
% disp(llik1);
bf=llik1-lmlik2;
disp(bf);
%bf=exp(llik1)/mlik2;
% else
%    bf=inf;
% end   