%Returns the wavelet reconstruction evaluated on the grid
%[0:2^-L:1-2^-L]^2, for a given vector of wavelet coefficients
%
%alpha: wavelet coefficients in a 2 dimensional grid.
%a: degree of Daubechies wavelets.
%L : resolution of reconstruction (up to 2^L).
%j0 : dilation factor of the given wavelet coefficients.
%p: number of iterations to compute boundary scaling fns



function R = get2DReconstruction(alpha, a, L, j0, p)
    sc = IWT2_CDJV_noP(real(alpha), ceil(log2(2*a))+1,a) + ...
            1i*IWT2_CDJV_noP(imag(alpha), ceil(log2(2*a))+1,a);
        
    if nargin < 5
        M = getScalingMatrix(a,j0, L);
    else
    	M = getScalingMatrix(a,j0, L, p);
    end
        
    R_temp = zeros(2^j0,2^L);
    for x=1:2^j0
        R_temp(x,:) = M*sc(x,:).';
    end
    R= zeros(2^L,2^L);
    for y = 1:2^L
        R(y,:) = M*R_temp(:,y);
    end
    
end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
