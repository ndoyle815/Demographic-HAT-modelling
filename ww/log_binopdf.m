%Computes the log beta-binomial pmf using N samples, with X data and with a
%probability distribution with a mean of p, and shape a (corresponds to
%parameter a in beta distribution)

function [LB]=log_binopdf(X,N,p)
if isnan(X) == 1
    LB =0;
elseif (N==0 && X==0) || (X==0 && p==0)
    LB=0; %probability of 1 of choosing 0 people from 0 people (LL of 0, exp(0)=1) or 0 people given zero probability

elseif (X~=0 && p==0) || (X~=0 && N==0)
            LB=-Inf; %impossible to choose people if probability is 0 or nobody is chosen

else 
        
        ix2=find(X>N);
        ix=find(X<=N);
        LB(ix)=gammaln(N+1) - gammaln(X(ix)+1) - gammaln(N-X(ix)+1) + log(p)*X(ix) + (N-X(ix))*log(1-p);
        
        LB(ix2)=-Inf*ones(1,length(ix2));
    
end
end