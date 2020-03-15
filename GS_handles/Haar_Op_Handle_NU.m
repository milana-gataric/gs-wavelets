% Applies 1D GS operator with Haar wavelets to a vector when samples are 
% taken nonuniformly in a form of a handle function y = G(wc,mode) ready to 
% be used in matlab's function 'lsqr' or 'spg_bp' from spgl1
%
% Output:  
%       y = G(wc,mode) is y = G *wc if mode == 1 (or 'notransp'), and
%                         y = G'*wc if mode == 2 (or 'transp'). 
%       where G is the generalized sampling operator 
%       G : wc [a vector of 2^R 1D Haar coefficients] -> y [a vector of M nonuniform Fourier coefficients]
% Input:
%       mode: argument of the function handle; 
%             mode == 1 (or 'notransp') computes the forward operation and 
%             mode == 2 (or 'transp') computes the adjoint operation
%       wc:   argument of the handle function;
%             a vector of 2^R Haar wavelet coefficients for mode == 1 and
%             a vector of M nonuniform Fourier coefficients for mode == 2
%       R: maximum scale of wavelet coefficients
%       ft_sca: Fourier transform of the scaling function
%       mu: a vector of M density compensation factors
%       st: an initial structure for NUFFT 
%       mask: an optional argument which represents a mask for the Fourier coefficients

function y = Haar_Op_Handle_NU(mode, wc, R, ft_sca, mu, st, mask)

L=2^R;

filter = [1 1] ./ sqrt(2);

if or(strcmp(mode,'notransp'),mode==1)
    x = IWT_PO(real(wc),0,filter)+1i*IWT_PO(imag(wc),0,filter);
    
    x_pad = zeros(L, 1);
    x_pad(1:2^R) = x;
    
    y = nufft(x_pad,st);

    y = ft_sca.*y;
    y = y.*sqrt(1/L);
    
    y = sqrt(mu).*y;
    
    if exist('mask','var')
        y = y.*mask;
    end
else
    if exist('mask','var')
        wc = wc.*mask;
    end
    wc = wc.*sqrt(1/L);
    wc = sqrt(mu).*wc;
    x=conj(ft_sca).*wc;
    
    x_fft = nufft_adj(x,st);
    
    x_trunc=x_fft(1:2^R);
    
    y =FWT_PO(real(x_trunc), 0,filter)+1i*FWT_PO(imag(x_trunc), 0,filter);


end

end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
