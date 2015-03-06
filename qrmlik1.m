function eul=qrmlik1(g0,g1,y,fy,alp,n,xm)

% b0=zeros(size(g0)); 
ul=zeros(length(g1),1);
% qy=zeros(length(y),length(g0)); u=zeros(length(y),length(g0));

b0=g0./(1-g0.^2);
b1=g1./(1-g1.^2);

%% For general double quadrature
% for i=1:length(g0)
%     for j=1:length(g1)
%         %     b0(i)=g0(i)/(1-g0(i)^2);
%         qy = b0(i) + b1(j)*fy';
%         u=y'-qy;
%         ul(i,j)= -n*log(sum(u.*(alp-(u<0))))-xm + gammaln(n) +...
%             n*log(alp*(1-alp))+log(1+g0(i)^2)+log(1+g1(j)^2)-2*log(1-g0(i)^2)-...
%             2*log(1-g1(j)^2);
%         %is there supposed to be a b1(i) over here?
%         ul(i,j)= ul(i) -log(2*pi) - 0.5*log(100^2) -0.5*sum([b0(i) b1(j)].^2)/100;
%     end
% end

%% For Matlab dblquad
for i=1:length(g0)
        %     b0(i)=g0(i)/(1-g0(i)^2);
        qy = b0(i) + b1(i)*fy';
        u=y'-qy;
        ul(i)= -n*log(sum(u.*(alp-(u<0))))-xm + gammaln(n) +...
            n*log(alp*(1-alp))-0.5*sum([b0(i) b1(i)].^2)/100;        
end
eul=ul+log(1+g0.^2)+log(1+g1.^2)-2*log(1-g0.^2)-2*log(1-g1.^2)...
                                               -log(2*pi) - 0.5*log(100^2);


% eul = exp(ul);
  

