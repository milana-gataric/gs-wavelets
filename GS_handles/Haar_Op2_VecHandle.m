% Applies 2D GS operator with Haar wavelets to a vector when samples 
% are taken uniformly in a form of a handle function h = G(mode,x) ready to 
% be used in matlab's function 'lsqr' or 'spg_bp' from spgl1
%
% Output:  
%       h = G(mode,x) is h = G *x if mode == 1 (or 'notransp');
%                        h = G'*x if mode == 2 (or 'transp')
%                        where G is the generalized sampling operator which
%                        maps a vector of (2^R)^2 2D Haar coefficients 
%                        to the vector of (2*S*2^R/eps)^2 
%                        2D UNIFORM Fourier coefficients
% Input:
%       mode:   argument of the function handle; 
%               mode == 1 (or 'notransp') computes the forward operation and 
%               mode == 2 (or 'transp') computes the adjoint operation
%       x:      argument of the handle function;
%               a vector of (2^R)^2 2D wavelet coefficients for mode == 1 and
%               a vector of (2*S*2^R/eps)^2 2D Fourier coefficients for mode == 2
%       S:      natural number, this determines the Fourier sampling range: -S*2^R/eps, ... , 2^R*S/eps-1
%       eps:    Fourier sampling density, require 1/eps to be a natural number
%       R:      maximum scale of wavelet coefficients
%       ft_sca: Fourier transform of the main, left, and right scaling function
%       mask:   an optional argument which represents a mask for the Fourier coefficients

function z = Haar_Op2_VecHandle(mode, x, S, eps, R, ft_sca, mask)
    
    N=sqrt(length(x)); 
    wc=reshape(x,N,N);    
    
    if nargin == 6
        z = Haar_Op2_Wavelet(mode, wc, S, eps, R, ft_sca);
    else
        z = Haar_Op2_Wavelet(mode, wc, S, eps, R, ft_sca, mask);
    end
    
    z=z(:);
end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
