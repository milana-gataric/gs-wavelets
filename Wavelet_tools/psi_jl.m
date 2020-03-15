%RETURNS: vector t, containing values of f(z) = sqrt(2^j)*psi(2^j*z -l) for
%z in S, and elements in S are uniformly spaced at 2^(-p).
function t=psi_jl(S,psi,p,j,l)
    n=ceil(2^(j+p)*S-l*2^(p))-2^p*min(S);
    n=max(1,min(n,length(psi))); %indices <1 will be set to 1, indices >length(psi) are set to length(psi)

        t=sqrt(2^j)*psi(n);
end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
