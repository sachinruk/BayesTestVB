function ul=qrmlik1_s(g0,g1,y,fy,alp,n)

% b0=zeros(size(g0)); 

% qy=zeros(length(y),length(g0)); u=zeros(length(y),length(g0));
% s=size(g0); g0=g0(:); g1=g1(:);
b0=g0./(1-g0.^2);
b1=g1./(1-g1.^2);
% 
% ul=zeros(length(g0),1);
% %% For general double quadrature
% for i=1:length(g0)
% %     qy = b0(i) + b1(i)*fy';
% %     u=y'-qy;
%     u=y-b0(i) - b1(i)*fy;
%     ul(i)= -n*log(sum(u.*(alp-(u<0))));        
% end
% 
% ul=ul+log(1+g0.^2)+log(1+g1.^2)-2*log(1-g0.^2)-2*log(1-g1.^2)...
%     -log(2*pi)-0.5*log(100^2)+ gammaln(n)-0.5*sum([b0 b1].^2,2)/100+...
%     n*log(alp*(1-alp));


% X=[ones(length(fy),1) fy'];
block_size=20; blocks=length(g0)/block_size; start=0;
b0_=reshape(b0,blocks,block_size); b1_=reshape(b1,blocks,block_size);
ul=zeros(1,length(g0));
for i=1:blocks
%     qy=bsxfun(@plus,fy*b1((start+1):(start+block_size))',...
%                                         b0((start+1):(start+block_size))');
%     u=bsxfun(@minus,y,qy);
    u=bsxfun(@minus,y,bsxfun(@plus,fy*b1_(i,:),b0_(i,:)));
    ul((start+1):(start+block_size))=sum(u.*(alp-(u<0)));
    start=start+block_size;
end

ul=-n*log(ul)+log(1+g0.^2)+log(1+g1.^2)-2*log(1-g0.^2)-2*log(1-g1.^2)...
    -log(2*pi)-0.5*log(100^2)+ gammaln(n)-0.5*sum([b0 b1].^2,2)/100+...
    n*log(alp*(1-alp));


% ul=reshape(ul,s);
% eul = exp(ul);