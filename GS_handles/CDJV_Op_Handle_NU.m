% Applies 1D GS operator with boundary wavelets to a vector in the 
% nonuniform case in a form of a handle function y = G(wc,mode) ready to be 
% used in matlab's function 'lsqr' or 'spg_bp' spgl1
%
% Output:  
%       y = G(wc,mode) is   y = G *wc          if mode == 1 (or 'notransp');
%                           y = adjoint(G)*wc  if mode == 2 (or 'transp')
%       where G is the generalized sampling operator 
%       G : wc [a vector of 2^R 1D wavelet coef] -> y [a vector of M nonuniform Fourier coef] 
%                         
%
% Input:
%       mode: argument of the function handle; 
%             mode == 1 (or 'notransp') computes the forward operation and 
%             mode == 2 (or 'transp') computes the adjoint operation
%       wc:   argument of the handle function;
%             a vector of 2^R 1D wavelet coefficients with a vanishing moments for mode == 1 and
%             a vector of M 1D Fourier coefficients for mode == 2
%       R: maximum scale of wavelet coefficients
%       a: number of vanishing moments
%       ft_sca, ft_sca_L, ft_sca_R: Fourier transform of the main, left, and right scaling function
%       sp: a vector of M nonuniform samples
%       mu: a vector of M density compensation factors
%       st: an initial structure for NUFFT
%       mask: an optional argument which represents a mask for the Fourier coefficients

function y = CDJV_Op_Handle_NU(mode, wc, R, a, ft_sca, ft_sca_L, ft_sca_R, sp, mu, st, mask)


if or(strcmp(mode,'notransp'),mode==1)
    
    x = IWT_CDJV_noP(real(wc),ceil(log2(2*a))+1,a) + 1i*IWT_CDJV_noP(imag(wc),ceil(log2(2*a))+1,a);
    if nargin==11
        y = CDJV_Op_Scaling_NU(mode, x, R, a, ft_sca, ft_sca_L, ft_sca_R, sp, mu, st, mask);
    else
        y = CDJV_Op_Scaling_NU(mode, x, R, a, ft_sca, ft_sca_L, ft_sca_R, sp, mu, st);
    end    

else
    
    if nargin==11
        sc = CDJV_Op_Scaling_NU(mode, wc, R, a, ft_sca, ft_sca_L, ft_sca_R, sp, mu, st, mask);
    else
        sc = CDJV_Op_Scaling_NU(mode, wc, R, a, ft_sca, ft_sca_L, ft_sca_R, sp, mu, st); 
    end
    y = FWT_CDJV_noP(real(sc), ceil(log2(2*a))+1,a)+1i*FWT_CDJV_noP(imag(sc), ceil(log2(2*a))+1,a);
    
end

end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
