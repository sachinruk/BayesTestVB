function eul=dq2mlik2(g0,g1,g2,y,fy,n,xm)

  for i=1:length(g0)
      b0(i)=g0(i)/(1-g0(i)^2);
  end    
  b1=g1/(1-g1^2);b2=g2/(1-g2^2);
  
  for i=1:length(g0)
      hit(:,i)=(y<fy)';
      hn=hit(2:n,i);hn1=hit(1:n-1,i);
      hitp(:,i) = b0(i) + b1*hn1 + b2*fy(2:n)';
      ul(i)= sum(hn.*hitp(:,i)) - sum(log(1+exp(hitp(:,i))))-xm+log(1+g0(i)^2)+log(1+g1^2)+log(1+g2^2)-2*log(1-g0(i)^2)-2*log(1-g1^2)-2*log(1-g2^2);
      ul(i)= ul(i) -3/2*log(2*pi) - 0.5*log(100^3) -0.5*[b0(i) b1 b2]*[b0(i);b1;b2]/100;
  end
  
  eul = exp(ul);