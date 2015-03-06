function eul=qrmlik2(g0,y,fy,alp,n)

  b0=g0(1)/(1-g0(1)^2);b1=g0(2)/(1-g0(2)^2);
  
  qy = b0 + b1*fy;
  u=y-qy;

  eul= -n*log(sum(u.*(alp-(u<0)))) +gammaln(n)+n*log(alp*(1-alp))+...
      log(1+g0(1)^2)+log(1+g0(2)^2)-2*log(1-g0(1)^2)-2*log(1-g0(2)^2);
  eul= -eul;
  %eul = exp(ul);