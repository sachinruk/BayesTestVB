function eul=dq2mlik2_s(g0,g1,g2,y,fy,n,xm)
s=size(g0); g0=g0(:); g1=g1(:); g2=g2(:);
b0=g0./(1-g0.^2);
b1=g1./(1-g1.^2);
b2=g2./(1-g2.^2);

hit=(y<fy)'; hn=hit(2:n); hn1=hit(1:n-1); fy_tail=fy(2:n)';

% c=b1*hn1 + b2*fy_tail;
% hitp=bsxfun(@plus,c(:,ones(length(b0),1)),b0);
% ul=sum(bsxfun(@times,hn,hitp)-log(1+exp(hitp)));
ul=zeros(size(g0));
for i=1:length(g0)
  hitp = b0(i) + b1(i)*hn1 + b2(i)*fy_tail;
  ul(i)= sum(hn.*hitp-log(1+exp(hitp)));
end
ul=ul+log(1+g0.^2)+log(1+g1.^2)+log(1+g2.^2)-2*log(1-g0.^2)...
    -2*log(1-g1.^2)-2*log(1-g2.^2)-3/2*log(2*pi) - 0.5*log(100^3) ...
    -0.5*sum([b0; b1; b2].^2)/100-xm;
eul = reshape(exp(ul),s);