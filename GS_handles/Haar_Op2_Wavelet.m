% Applies 2D GS operator with Haar wavelets to a matrix when samples 
% are taken uniformly in a form of a handle function y = G(mode,wc) 
%
% Output:  
%                G : wc [2^R times 2^R] -> y [2M times 2M] if mode == 1 (or 'notransp');
%       adjoint(G) : wc [2M times 2M] -> y [2^R times 2^R] if mode == 2 (or 'transp'). 
% Input:
%       mode: mode == 1 (or 'notransp') computes the forward operation and 
%             mode == 2 (or 'transp') computes the adjoint operation
%       wc:   a matrix 2^R times 2^R of 2D HAAR coefficients for mode == 1 and
%             a matrix 2M times 2M (M=S*2^R/eps) of 2D UNIFORM Fourier coefficients for mode == 2
%       S: natural number, this determines the Fourier sampling range: -S*2^R/eps, ... , 2^R*S/eps-1
%       eps: Fourier sampling density, require 1/eps to be a natural number
%       R: maximum scale of wavelet coefficients
%       ft_sca: Fourier transform of the main, left, and right scaling function
%       mask: an optional argument which represents a mask for the Fourier coefficients


function y_concat = Haar_Op2_Wavelet(mode, wc, S, eps, R, ft_sca, mask)
    
    L=2^R/eps;
    filter = [1 1] ./ sqrt(2);
    
    if or(strcmp(mode,'notransp'),mode==1)
        
        x = IWT2_PO(real(wc),0,filter)+1i*IWT2_PO(imag(wc),0,filter);
        y_concat = zeros(2*S*L);

        for iy = 1:2^R
            y_concat(:,iy) = Haar_Op_Handle_Scaling(mode, x(:,iy), S, eps, R, ft_sca);
        end
        for ix = 1:2*S*L
            y_concat(ix,:) =  Haar_Op_Handle_Scaling(mode, y_concat(ix,1:2^R).', S, eps, R, ft_sca);
        end
        
        if nargin == 7
            y_concat = mask.*y_concat;
        end
        
    else
        if nargin == 7
            wc = mask.*wc;
        end
        
        y_concat = zeros(2^R);
        y_interm = zeros(2*S*L, 2^R);
        
        
        for ix = 1:2*S*L
            y_interm(ix,:) =  Haar_Op_Handle_Scaling(mode, wc(ix,:).', S, eps, R, ft_sca);
        end
        
        for iy = 1:2^R
            y_concat(:,iy) = Haar_Op_Handle_Scaling(mode, y_interm(:,iy), S, eps, R, ft_sca);
        end
        x_trunc = y_concat;
        y_concat = FWT2_PO(real(x_trunc),0,filter)+1i*FWT2_PO(imag(x_trunc),0,filter);
        
    end
end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
