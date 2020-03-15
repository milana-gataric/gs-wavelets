% This function is the same as function 'CDJV_Op_Handle', except that now 
% wc represents 1D boundary corrected scaling coefficients instead of 1D 
% boundary corrected wavelet coefficients

function y_concat = CDJV_Op_Scaling(mode, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask)

L=2^R/eps;
ft_sca_R = (ft_sca_R.* repmat(exp(2*pi*1i *eps *(-S*L:S*L-1)'),1,a));

if or(strcmp(mode,'notransp'),mode==1)
    
    x_pad = zeros(L, 1);
    x_pad(2:2^R-2*a+1) = wc(a+1:2^R-a);
    y= L*ifft(x_pad);
    
    if S<1
        y = fftshift(y);
    end
    
    y_concat = repmat(y,2*S,1);

    y_concat = ft_sca.*y_concat;


    y_concat = y_concat + ft_sca_L*wc(1:a) + ft_sca_R * wc(2^R-a+1:2^R);


    y_concat = y_concat.*sqrt(1/L);
    if nargin == 10
        y_concat = mask.*y_concat;
    end

else
    if nargin == 10
        wc = mask.*wc;
    end
    wc=L^(-1/2)*wc;
    x=conj(ft_sca).*wc;
    
    if S<1
        x = fftshift(x);
    end
    
    x_fft=zeros(L,1);
    for j=1:2*S
        x_fft=x_fft+ fft(x(j*L-L+1:j*L));
    end
    x_trunc=zeros(2^R,1);
    
    x_trunc(a+1:2^R-a) = x_fft(2:2^R-2*a+1);
    
    x_trunc(1:a) = ft_sca_L'*wc;
    x_trunc(2^R-a+1:2^R) = ft_sca_R'*wc;
    
    
    y_concat = x_trunc;
    
end

end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
