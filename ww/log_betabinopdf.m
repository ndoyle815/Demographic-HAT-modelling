%Computes the log beta-binomial pmf using N samples, with X data and with a
%probability distribution with a mean of pi, and correlation coeffecient rho 

function [LB]=log_betabinopdf(X,N,p,rho)

if isnan(X) == 1
    LB =0;
elseif (N==0 && X==0) || (X==0 && p==0)
    LB=0; %probability of 1 of choosing 0 people from 0 people (LL of 0, exp(0)=1) 

elseif (X~=0 && p==0) || (X~=0 && N==0)
            LB=-Inf; %impossible to choose people if probability is 0 or nobody is chosen

elseif rho==0 %no variance in probability => binomial distribution
    
        ix2=find(X>N);
        ix=find(X<=N);
        LB(ix)=gammaln(N+1) - gammaln(X(ix)+1) - gammaln(N-X(ix)+1) + log(p)*X(ix) + (N-X(ix))*log(1-p);
        
        LB(ix2)=-Inf*ones(1,length(ix2));
    
else 
        %computes b to give beta distribution with mean p given specific a
        a = p*(1/rho - 1);
        b=a*(1-p)/p;

        ix2=find(X>N);
        ix=find(X<=N);
        LB(ix)=gammaln(N+1)+gammaln(X(ix)+a) +gammaln(N-X(ix)+b) +gammaln(a+b)-gammaln(X(ix)+1)-gammaln(N-X(ix)+1)-gammaln(N+a+b)-gammaln(a)-gammaln(b);
        LB(ix2)=-Inf*ones(1,length(ix2));
    
end

end