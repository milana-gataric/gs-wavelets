% Applies 1D GS operator with boundary wavelets to a vector in a form of a 
% handle function y = G(wc,mode) ready to be used in matlab's function
% 'lsqr' or 'spg_bp' from spgl1 package
%
% Output:  
%       y = G(wc,mode) is   y = G *wc          if mode == 1 (or 'notransp');
%                           y = adjoint(G)*wc  if mode == 2 (or 'transp')
%       where G is the generalized sampling operator 
%       G : wc [a vector of 2^R 1D wavelet coef] -> y [a vector of 2*S*2^R/eps 1D uniform Fourier coef] 
%                         
%
% Input:
%       mode: argument of the function handle; 
%             mode == 1 (or 'notransp') computes the forward operation and 
%             mode == 2 (or 'transp') computes the adjoint operation
%       wc:   argument of the handle function;
%             a vector of 2^R 1D wavelet coefficients with a vanishing moments for mode == 1 and
%             a vector of 2*S*2^R/eps 1D Fourier coefficients for mode == 2
%       S: natural number, this determines the Fourier sampling range: -S*2^R/eps, ... , 2^R*S/eps-1
%       eps: Fourier sampling density, require 1/eps to be a natural number
%       R: maximum scale of wavelet coefficients
%       a: number of vanishing moments
%       ft_sca, ft_sca_L, ft_sca_R: Fourier transform of the main, left, and right scaling function
%       mask: an optional argument which represents a mask for the Fourier coefficients

function y = CDJV_Op_Handle(mode, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask)

if or(strcmp(mode,'notransp'),mode==1)
    
    x = IWT_CDJV_noP(real(wc),ceil(log2(2*a))+1,a) + 1i*IWT_CDJV_noP(imag(wc),ceil(log2(2*a))+1,a);
    if nargin == 10
        y = CDJV_Op_Scaling(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
    else
        y = CDJV_Op_Scaling(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
    end

else
    if nargin == 10
        sc = CDJV_Op_Scaling(mode, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
    else
        sc = CDJV_Op_Scaling(mode, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
    end
    y = FWT_CDJV_noP(real(sc), ceil(log2(2*a))+1,a)+1i*FWT_CDJV_noP(imag(sc), ceil(log2(2*a))+1,a);
    
end

end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
