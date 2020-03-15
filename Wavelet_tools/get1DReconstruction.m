% Returns the wavelet reconstruction evaluated at the
% points 0:2^-L:1-2^-L, for a given vector of wavelet coefficients
%
% alpha: wavelet coefficients, should be a column vector
%     a: degree of Daubechies wavelets
%     L: resolution of reconstruction (up to 2^L)
%    j0: dilation factor of the given wavelet coefficients
%     p: number of iterations to compute boundary scaling fns

function R = get1DReconstruction(alpha, a, L, j0, p)
	sc = IWT_CDJV_noP(real(alpha), ceil(log2(2*a))+1,a) + ...
            1i*IWT_CDJV_noP(imag(alpha), ceil(log2(2*a))+1,a);
    if nargin < 5    
        M = getScalingMatrix(a,j0, L);
    else
        M = getScalingMatrix(a,j0, L, p);
    end
              
    R = M*sc;
    
    
end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
