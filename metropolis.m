function [beta, AR]=metropolis(y,X,alpha,stepSize)
samples=10000;
%intial beta sample
%let variance be 1/10th of q(beta)
beta=zeros(2,5000);
beta(:,1)=randn(2,1);
l_pby_old=pby(y,X,beta(:,1),alpha);
accept=0;
for i=1:samples
    beta_new=transition_q(beta(:,i),stepSize);
    l_pby=pby(y,X,beta_new,alpha);
    if l_pby>l_pby_old
        beta(:,i+1)=beta_new;
        l_pby_old=l_pby;
        accept=accept+1;
    else
        if exp(l_pby-l_pby_old)>rand
            beta(:,i+1)=beta_new;
            l_pby_old=l_pby;
            accept=accept+1;
        else
            beta(:,i+1)=beta(:,i);
        end
    end
end

AR=accept/samples;
burn_in=2000;
beta=beta(:,(burn_in+1):end)';

% %log likelihood
%  lik=zeros(samples-burn_in,1); 
% for i=(burn_in+1):samples
%     lik(i)=pby(y,X,beta(:,i),alpha);
% end
% llik=logsumexp(lik)-log(samples-burn_in);
    

function l_pby=pby(y,X,beta,alpha)

n=length(y); d=length(beta);
error=y-X*beta;
I=error<0;
l_pb_y=n*log(alpha)+n*log(1-alpha)+gammaln(n)-n*log(sum(error.*(alpha-I)));
l_pb=-d/2*log(2*pi*1000)-1/(2*1000)*(beta'*beta);

l_pby=l_pb_y+l_pb;

function beta=transition_q(beta_t_1,stepSize)
d=length(beta_t_1);
beta=beta_t_1+stepSize*randn(d,1);