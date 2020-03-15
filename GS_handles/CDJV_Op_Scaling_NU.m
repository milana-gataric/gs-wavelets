% This function is the same as function 'CDJV_Op_Handle_NU', except that now 
% wc represents 1D boundary corrected scaling coefficients instead of 1D 
% boundary corrected wavelet coefficients

function y = CDJV_Op_Scaling_NU(mode, wc, R, a, ft_sca, ft_sca_L, ft_sca_R, SamplePoints, mu, st, mask)

L=2^R;
ft_sca_R = (ft_sca_R.* repmat(exp(-2*pi*1i *SamplePoints),1,a));


if or(strcmp(mode,'notransp'),mode==1)
    
    x_pad = zeros(2^R, 1);
    x_pad(2:2^R-2*a+1) = wc(a+1:2^R-a);
    
    y_concat = nufft(x_pad,st);

    y_concat = ft_sca.*y_concat;

    y_concat = y_concat + ft_sca_L*wc(1:a) + ft_sca_R * wc(2^R-a+1:2^R);

    y_concat = y_concat.*sqrt(1/L);

    y = sqrt(mu).*y_concat;
    
    if nargin == 11
        y = mask.*y;
    end
    
else
    wc = L^(-1/2)*wc;
    wc = sqrt(mu).*wc;
    if nargin == 11
        wc = mask.*wc;
    end
    
    x = conj(ft_sca).*wc;
    
    x_fft = nufft_adj(x,st);
    
    x_trunc = zeros(2^R,1);
    
    x_trunc(a+1:2^R-a) = x_fft(2:2^R-2*a+1);
    
    x_trunc(1:a) = ft_sca_L'*wc;
    x_trunc(2^R-a+1:2^R) = ft_sca_R'*wc;
    
    
    y = x_trunc;
    
end

end


% Copyright (c) 2014. Clarice Poon and Milana Gataric
