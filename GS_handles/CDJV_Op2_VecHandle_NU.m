% Applies 2D GS operator with boundary wavelets to a vector when samples 
% are taken nonuniformly a form of a handle function z = G(mode,x) ready 
% to be used in matlab's function 'lsqr' or 'spg_bp' from spgl1
%
% Output:  
%       z = G(mode,x) is  z = G*x           if mode == 1 (or 'notransp');
%                         z = adjoint(G)*x  if mode == 2 (or 'transp')
%       where G is the generalized sampling operator 
%       G : x [a vector of (2^R)^2 2D wavelet coeff] -> y [a vector of M nonuniform Fourier samples] 
%
% Input:
%       mode: argument of the function handle; 
%             mode == 1 (or 'notransp') computes the forward operation and 
%             mode == 2 (or 'transp') computes the adjoint operation
%       x:    argument of the handle function;
%             a vector of (2^R)^2 2D wavelet coefficients for mode == 1 and
%             a vector of length(sp) NONUNIFORM Fourier coefficients for mode == 2
%       R: maximum scale of wavelet coefficients
%       a: number of vanishing moments
%       ft_sca_x, ft_sca_L_x, ft_sca_R_x: Fourier transform of the main, left, 
%                                         and right scaling function at sp(:,1)
%       ft_sca_y, ft_sca_L_y, ft_sca_R_y: Fourier transform of the main, left, 
%                                         and right scaling function at sp(:,2)
%       sp: a M times 2 matrix of sampling points
%       mu: a vector of M density compensetaion factors
%       st, stx, sty: initialize structure for NUFFT
%       mask: an optional argument which represents a mask for the Fourier coefficients

function z = CDJV_Op2_VecHandle_NU(x, mode, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty, mask )
    if or(strcmp(mode,'notransp'),mode==1)
        N=sqrt(length(x));    
        wc=reshape(x,N,N); 
        if nargin == 15
            z = CDJV_Op2_Handle_NU(mode, wc, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
        else
            z = CDJV_Op2_Handle_NU(mode, wc, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty, mask);
        end
    else
        if nargin == 15
            z = CDJV_Op2_Handle_NU(mode, x, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
        else
            z = CDJV_Op2_Handle_NU(mode, x, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty, mask);
        end
        z = z(:);
    end
end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
