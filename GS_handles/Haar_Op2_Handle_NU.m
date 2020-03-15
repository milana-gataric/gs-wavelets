% Applies 2D GS operator with Haar wavelets to a matrix when samples 
% are taken nonuniformly in a form of a handle function y = G(mode,wc) 
%
% Output:  
%       y = G(mode,x) is  y = G*wc           if mode == 1 (or 'notransp');
%                         y = adjoint(G)*wc  if mode == 2 (or 'transp')
%       where G is the generalized sampling operator 
%       G : wc [a matrix of 2^R times 2^R of 2D Haar coeff] -> y [a vector of M nonuniform Fourier samples] 
%
% Input:
%       mode: argument of the function handle; 
%             mode == 1 (or 'notransp') computes the forward operation and 
%             mode == 2 (or 'transp') computes the adjoint operation
%       x:    argument of the handle function;
%             a vector of (2^R)^2 2D HAAR coefficients for mode == 1 and
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

function y = Haar_Op2_Handle_NU(mode, wc, R, ft_sca_x, ft_sca_y, mu, st, mask)

    L=2^R;
    filter = [1 1] ./ sqrt(2);
    
    if or(strcmp(mode,'notransp'),mode==1)
       
        x = IWT2_PO(real(wc),0,filter)+1i*IWT2_PO(imag(wc),0,filter);
        xt = x.';
        
        y = nufft(xt,st);

        
        y = ft_sca_x.*(ft_sca_y.*y);
        
        y = y.*sqrt(mu)*1/L;
        
        if nargin == 8
            y = mask.*y;
        end

    else

        wc=wc./L;
        wc = sqrt(mu).*wc;
        if nargin == 8
            wc = mask.*wc;
        end
        
        x = conj(ft_sca_y).*(conj(ft_sca_x).*wc);        
        x_trunc = nufft_adj(x,st);
        x_trunc = x_trunc.';
       
        y = FWT2_PO(real(x_trunc),0,filter)+1i*FWT2_PO(imag(x_trunc),0,filter);

    end


end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
