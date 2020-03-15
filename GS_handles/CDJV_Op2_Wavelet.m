% Applies 2D GS operator with boundary wavelets to a matrix when samples 
% are taken uniformly in a form of a handle function y = G(mode,wc) 
%
% Output:  
%                G*wc : wc [a matrix 2^R times 2^R] -> y [a matrix 2M times 2M] if mode == 1 (or 'notransp');
%       adjoint(G)*wc : wc [a matrix 2M times 2M] -> y [a matrix 2^R times 2^R] if mode == 2 (or 'transp'). 
% Input:
%       mode: mode == 1 (or 'notransp') computes the forward operation and 
%             mode == 2 (or 'transp') computes the adjoint operation
%       wc:   a matrix 2^R times 2^R of 2D wavelet coefficients for mode == 1 and
%             a matrix 2M times 2M (M=S*2^R/eps) of 2D UNIFORM Fourier coefficients for mode == 2
%       S: natural number, this determines the Fourier sampling range: -S*2^R/eps, ... , 2^R*S/eps-1
%       eps: Fourier sampling density, require 1/eps to be a natural number
%       R: maximum scale of wavelet coefficients
%       a: number of vanishing moments
%       ft_sca, ft_sca_L, ft_sca_R: Fourier transform of the main, left, and right scaling function
%       mask: an optional argument which represents a mask for the Fourier coefficients


function y_concat = CDJV_Op2_Wavelet(mode, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask)
    L=2^R/eps;
    
    if or(strcmp(mode,'notransp'),mode==1)
        x = IWT2_CDJV_noP(real(wc),ceil(log2(2*a))+1,a) +1i * IWT2_CDJV_noP(imag(wc),ceil(log2(2*a))+1,a);
        y_concat = zeros(2*S*L);

        for iy = 1:2^R
            y_concat(:,iy) = CDJV_Op_Scaling(mode, x(:,iy), S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
        end

        for ix = 1:2*S*L
            y_concat(ix,:) =  CDJV_Op_Scaling(mode, y_concat(ix,1:2^R).', S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
        end
        
        if nargin == 10
            y_concat = mask.*y_concat;
        end
        
    else
        if nargin == 10
            wc = mask.*wc;
        end
        
        y_concat = zeros(2^R);
        y_interm = zeros(2*S*L, 2^R);
        
        
        for ix = 1:2*S*L
            y_interm(ix,:) =  CDJV_Op_Scaling(mode, wc(ix,:).', S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
        end
        
        for iy = 1:2^R
            y_concat(:,iy) = CDJV_Op_Scaling(mode, y_interm(:,iy), S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
        end
        x_trunc = y_concat;
        y_concat = FWT2_CDJV_noP(real(x_trunc), ceil(log2(2*a))+1,a)+1i*FWT2_CDJV_noP(imag(x_trunc), ceil(log2(2*a))+1,a);
        
    end
end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
