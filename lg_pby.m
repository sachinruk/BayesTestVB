function l_pby=lg_pby(beta,y,X,alpha)

n=length(y); d=length(beta);
error=y-X*beta;
I=error<0;
l_pb_y=n*log(alpha)+n*log(1-alpha)+gammaln(n)-n*log(sum(error.*(alpha-I)));
l_pb=-d/2*log(2*pi*1000)-1/(2*1000)*(beta'*beta);

l_pby=l_pb_y+l_pb;