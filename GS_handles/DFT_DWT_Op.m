% Description:
% Function handle to test compressed sensing using DFT measurements and
% using DWT as the sparsifying transform.
%
% Inputs:
% mode: set to 1 for forward operation, anything else for adjoint operation
% wc: input vector of length 2^R
% a: number of vanishing moments of wavelet tranform
% R: resolution of input
%
%Output:
% y = U*wc or U'*wc, where U= DFT DWT^{-1}
%
function y = DFT_DWT_Op(mode, wc, a, R)
N = sqrt(length(wc));
wc = reshape(wc, N,N);
if mode==1 %forward
    
    if a==1
        filter = [1 1] ./ sqrt(2);
        x = IWT2_PO(real(wc),0,filter)+1i*IWT2_PO(imag(wc),0,filter);
    else
        x = IWT2_CDJV(real(wc),ceil(log2(2*a))+1,a) + 1i*IWT2_CDJV(imag(wc),ceil(log2(2*a))+1,a);
    end
    
    y =  R.*ifft2(x)*N;
    y = y(:);
else %transpose
    sc = fft2(R.*wc)/N;
    
    if a==1
        filter = [1 1] ./ sqrt(2);
        y = FWT2_PO(real(sc), 0,filter)+1i*FWT2_PO(imag(sc), 0,filter);
    else
        y = FWT2_CDJV(real(sc), ceil(log2(2*a))+1,a)+1i*FWT2_CDJV(imag(sc), ceil(log2(2*a))+1,a);
    end
    y=y(:);
end

end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
