function eul=dq2mlik22b(g0,y,fy,n)
%dq2mlik22b(b2,y,fy,n)
  b0=g0(1);b1=g0(2);b2=g0(3);
  
  hit=(y<fy);
  
  hitp = b0 + b1*hit(1:n-1)+b2*fy(2:n);
  
  eul= sum(hit(2:n).*hitp)-sum(log(1+exp(hitp)));%+log(1+g0(1)^2)+log(1+g0(2)^2)+log(1+g0(3)^2)-2*log(1-g0(1)^2)-2*log(1-g0(2)^2)-2*log(1-g0(3)^2);
  eul=-eul;
  %eul = exp(ul);