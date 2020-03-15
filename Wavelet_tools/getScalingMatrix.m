%Returns a matrix containing pointwise evaluations (by rows) of the scaling
%functions (by columns) of dilation factor 2^j0, at the points 0:2^-L:1-2^-L.
%The first a columns of M are the left boundary scaling functions.
%The last a columns of M are the right boundary scaling functions.

%p is optional, the number of iterations to compute boundary scaling fns

function M = getScalingMatrix(a, j0, L, p)

    %evaluations of the internal scaling functions are calculated using 
    %the dyadic dilation algorithm, so should be exact.
    [S,Phi] =  getScalingFn(a, L);
    
    S = S-a+1;
    
    ind0 = find(S==0);
    trunc = @(A) A(ind0:ind0+2^L-1);
    
    fun = @(k) trunc(psi_jl(S,Phi,L,j0,k));
    
    t = a:2^j0-a-1;
    M = cell2mat(arrayfun(fun,t,'UniformOutput', false));

    %%evaluations of boundary scaling functions are calculated using the
    %%cascade algorithm, so non-exact.
    if nargin == 4
        n=2^p;
    else 
        L0 = max(L, 20);
        n=2^L0;
    end
    
    M_left = zeros(2^L, a);
    M_right = zeros(2^L, a);
    for l=a:-1:1
        wav = MakeWavelet(j0,l,'CDJV',a,'Father',n)*sqrt(n); 
        wav = wav(1:n/2^L:end);
        M_left(:,l) = wav.';
    end
    
    for l=2^j0-a+1:2^j0
        wav = MakeWavelet(j0,l,'CDJV',a,'Father',n)*sqrt(n); 
        wav = wav(1:n/2^L:end);
        M_right(:,l-2^j0+a) =wav.';
    end
    M = horzcat(M_left, M, M_right);

end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
