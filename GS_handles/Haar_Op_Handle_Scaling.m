% This function is the same as function 'Haar_Op_Handle', except that now 
% wc represents 1D scaling coefficients instead of 1D wavelet coefficients

function y_concat = Haar_Op_Handle_Scaling(mode, wc, S, eps, R, ft_sca, mask)

L=2^R/eps;


if or(strcmp(mode,'notransp'),mode==1)

    x_pad = zeros(L, 1);
    x_pad(1:2^R) = wc;
    y= L*ifft(x_pad);
    y_concat = repmat(y,2*S,1);

    y_concat = ft_sca.*y_concat;
    y_concat = y_concat.*sqrt(eps/2^R);
    
    if exist('mask', 'var')
        y_concat = y_concat.*mask;
    end
else
    if exist('mask', 'var')
        wc = wc.*mask;
    end
    wc = wc.*sqrt(1/L);
    x=conj(ft_sca).*wc;
    
    x_fft=zeros(L,1);
    for j=1:2*S
        x_fft=x_fft+ fft(x(j*L-L+1:j*L));
    end
    y_concat=x_fft(1:2^R);
    


end

end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
