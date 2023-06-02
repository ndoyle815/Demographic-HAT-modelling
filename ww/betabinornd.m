% drawn random beta binomial numbers

%drawn n samples with a
%probability distribution with a mean of p, and shape a (corresponds to
%parameter alpha in beta distribution)

function R= betabinornd(n,p,rho)%,M,N)

%If rho is a scalar
if length(rho)==1
    
    if rho==0
        R=binornd(n,p);
        
    else
        
        a = p.*(1./rho - 1);
        b=a.*(1-p)./p;
        
        R=binornd(n,betarnd(a,b));
        
    end

%If rho is a matrix of the same size as n and p
else
    
   %if overdispersion is zero for any entries then use binornd
    
    for ix=1:size(rho,1)
        
        iy=find(rho(ix,:)==0);
        R(ix,iy)=binornd(n(ix,iy),p(ix,iy));
    end
    
    %If overdispersion is non-zero then do the betabinornd calculation
    

    for ix2=1:size(rho,1)
        
        iy2=find(rho(ix2,:)~=0);
        a(ix2,iy2) = p(ix2,iy2).*(1./rho(ix2,iy2) - 1);
        b(ix2,iy2) = a(ix2,iy2).*(1-p(ix2,iy2))./p(ix2,iy2);

        R(ix2,iy2)=binornd(n(ix2,iy2),betarnd(a(ix2,iy2),b(ix2,iy2)));
    
    end

end