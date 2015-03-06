function  bf=bfdq2mlik2(y,fy,alp)

n=length(y);

% g0=-0.9:0.02:0;ng=length(g0);
% g1=g0;g2=g1;
% for i=1:ng
%     for j=1:ng
%        for k=1:ng 
%          lmlik(i,j,k)=dq2mlik(g0(i),g1(j),g2(k),y,fy,alp,n);
%        end
%     end   
% end

%mlik=exp(lmlik-max(max(lmlik)));

%marginal likelihood, double integral

%xm=max(max(max(lmlik)))
%[X,FVAL,EXITFLAG,OUTPUT] = fminunc(FUN,X0,...)
%A=[eye(3);-1*eye(3)];B=[0.99 0.99 0.99 0.99 0.99 0.99]';
% X0=[-0.2 0.2 -0.2]';
% [g,xm]=fmincon(@(g2) dq2mlik22(g2,y,fy,n),X0,A,B);
X0=[-1 -0.2 0.2];
[~,xm]=fminsearch(@(b2) dq2mlik22b(b2,y,fy,n),X0);%,A,B);
%marginal likelihood, double integral

xm=-xm;

%mlik2 = dblquad(@(x0,x1) dqmlik1(x0,x1,y,fy,alp,n,xm),-0.99,0.99,-0.99,0.99);
% mlik2 = triplequad(@(x0,x1,x2) dq2mlik2_s2(x0,x1,x2,y,fy,n,xm),-0.99,0.99,-0.99,0.99,-0.99,0.99);
% disp(mlik2);
mlik2 = triplequad(@(x0,x1,x2) dq2mlik2(x0,x1,x2,y,fy,n,xm),-0.99,0.99,-0.99,0.99,-0.99,0.99);
disp(mlik2);
% integral3(@(x0,x1,x2) dq2mlik2_s(x0,x1,x2,y,fy,n,xm),-0.99,0.99,-0.99,0.99,-0.99,0.99)
% dx=0.01; dy=dx; dz=dy; lim=0.999; lim_z=-0.999; range_z=0.05; k=1;
% a=zeros(10,1);
% while lim_z<lim
%     [X0, X1, X2]=meshgrid(-lim:dx:lim,-lim:dy:lim, ...
%                                           lim_z:dy:min(lim,lim_z+range_z));
%     a(k)=logsumexp(dq2mlik2_s(X0(:),X1(:),X2(:),y,fy,n,xm));    
%     disp(k);
%     k=k+1;
%     lim_z=lim_z+range_z+dz;    
%     disp(k);
% end
% mlik2_ =exp(logsumexp(a)+log(dx)+log(dy)+log(dz));
%mlik3 = quad(@(x2) dbint(x2,y,fy,alp,n,xm),-0.99,0.99);
%Bayes factor

%likelihood  b0=0, b1=0;

hit=(y<fy);
hitp = log(alp/(1-alp)).*ones(n-1,1);
llik = sum(hit(2:n).*hitp') - sum(log(1+exp(hitp)));

llik1=llik-xm;

bf=exp(llik1-log(mlik2));
%bf=exp(llik1)/mlik2;
