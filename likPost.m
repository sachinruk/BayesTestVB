function L=likPost(y,x,alpha,mu,Sigma)

n=10000;
%simulate n betas from mu,Sigma
L=chol(Sigma)';
beta=bsxfun(@plus,L*randn(size(x,2),n),mu);

l_lik=zeros(n,1);
for i=1:n
    l_pby=pby(y,x,beta(:,i),alpha);
    l_qby=qby(beta(:,i),mu,L);
    l_lik(i)=l_pby-l_qby;
end
L=logsumexp(l_lik)-log(n);


function l_pby=pby(y,X,beta,alpha)

n=length(y); d=length(beta);
error=y-X*beta;
I=error<0;
l_pb_y=n*log(alpha)+n*log(1-alpha)+gammaln(n)-n*log(sum(error.*(alpha-I)));
l_pb=-d/2*log(2*pi*1000)-1/(2*1000)*(beta'*beta);

l_pby=l_pb_y+l_pb;

function l_qby=qby(beta,mu,L)

e=beta-mu; err=L\e;
l_qby=-0.5*(log(2*pi)+err'*err)-sum(log(diag(L)));