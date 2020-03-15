% Returns the wavelet reconstruction evaluated on the grid
% 0:2^-p:1-2^-p, for a given vector of wavelet coefficients
%
% z: Haar wavelet coefficients.
% p : resolution of reconstruction (up to 2^p).
% R : dilation factor of the given wavelet coefficients.
function recons = get1DHaarReconstruction(z, p, R)
    filter = [1 1] ./ sqrt(2);
    sc = IWT_PO(real(z),0,filter)+1i* IWT_PO(imag(z),0,filter);
    f_recon = @(x) sc(floor(x*2^R)+1);
    recons = f_recon(0:2^-p:1-2^-p)*sqrt(2^R);
end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
