% Applies 1D GS operator with Haar wavelets to a vector when sampling 
% Fourier transform uniformly in a form of a handle function 
% y_concat = G(wc,mode) ready to be used in the matlab function 
%'lsqr' or 'spg_bp' from spgl1
%
% Output:  
%       y_concat = G(wc,mode) is y_concat = G *wc if mode == 1 (or 'notransp');
%                                y_concat = G'*wc if mode == 2 (or 'transp'). 
%                                where G is the generalized sampling operator which
%                                maps a vector of 2^R 1D Haar coefficients 
%                                to the vector of 2*S*2^R/eps 1D UNIFORM
%                                Fourier coefficients
% Input:
%       mode: argument of the function handle; 
%             mode == 1 (or 'notransp') computes the forward operation and 
%             mode == 2 (or 'transp') computes the adjoint operation
%       wc:   argument of the handle function;
%             a vector of Haar wavelet coefficients for mode == 1 and
%             a vector of Fourier coefficients for mode == 2
%       S: natural number, this determines the Fourier sampling range: -S*2^R/eps, ... , 2^R*S/eps-1
%       eps: Fourier sampling density, require 1/eps to be a natural number
%       R: maximum scale of wavelet coefficients
%       ft_sca: Fourier transform of the scaling function
%       mask: an optional argument which represents a mask for the Fourier coefficients

function y_concat = Haar_Op_Handle(mode, wc, S, eps, R, ft_sca, mask)

L=2^R/eps;

filter = [1 1] ./ sqrt(2);

if or(strcmp(mode,'notransp'),mode==1)
    x = IWT_PO(real(wc),0,filter)+1i*IWT_PO(imag(wc),0,filter);
    x_pad = zeros(L, 1);
    x_pad(1:2^R) = x;
    y= L*ifft(x_pad);
    
    if S<1
        y = fftshift(y);
    end
    
    y_concat = repmat(y,2*S,1);

    y_concat = ft_sca.*y_concat;
    y_concat = y_concat.*sqrt(eps/2^R);
    
    if exist('mask','var')
        y_concat = y_concat.*mask;
    end
else
    if exist('mask','var')
        wc = wc.*mask;
    end
    wc = wc.*sqrt(1/L);
    x=conj(ft_sca).*wc;
    
    if S<1
        x = fftshift(x);
    end
    
    x_fft=zeros(L,1);
    for j=1:2*S
        x_fft=x_fft+ fft(x(j*L-L+1:j*L));
    end
    
    
    
    x_trunc=x_fft(1:2^R);
    
    y_concat =FWT_PO(real(x_trunc), 0,filter)+1i*FWT_PO(imag(x_trunc), 0,filter);


end

end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
