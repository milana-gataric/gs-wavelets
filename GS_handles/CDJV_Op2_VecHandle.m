% Applies 2D GS operator with boundary wavelets to a vector when samples 
% are taken uniformly in a form of a handle function h = G(mode,x) ready to 
% be used in matlab's function 'lsqr' or 'spg_bp' from spgl1 package 
%
% Output:  
%       h = G(mode,x) is  h = G*x          if mode == 1 (or 'notransp');
%                         h = adjoint(G)*x if mode == 2 (or 'transp')
%       where G is the generalized sampling operator 
%       G : x [a vector of (2^R)^2 2D wavelet coeff] -> y [a vector of (2*S*2^R/eps)^2 uniform samples] 
%
% Input:
%       mode: argument of the function handle; 
%             mode == 1 computes the forward operation and 
%             mode == 2 computes the adjoint operation
%       x:    argument of the handle function;
%             a vector of (2^R)^2 2D wavelet coefficients for mode == 1 and
%             a vector of (2*S*2^R/eps)^2 2D UNIFORM Fourier coefficients for mode == 2
%       S: natural number, this determines the Fourier sampling range: -S*2^R/eps, ... , 2^R*S/eps-1
%       eps: Fourier sampling density, require 1/eps to be a natural number
%       R: maximum scale of wavelet coefficients
%       a: number of vanishing moments
%       ft_sca, ft_sca_L, ft_sca_R: Fourier transform of the main, left, and right scaling function
%       mask: an optional argument which represents a mask for the Fourier coefficients


function h = CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask)
    
    N=sqrt(length(x)); 
    wc=reshape(x,N,N); 
    
    if nargin == 9
        h = CDJV_Op2_Wavelet(mode, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
    else
        h = CDJV_Op2_Wavelet(mode, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
    end
    
    h=h(:);

end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
