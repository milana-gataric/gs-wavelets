%%Returns Haar wavelet reconstruction evaluated at the
%%points 0:2^-p:1-2^-p, for a given vector of wavelet coefficients
%
%wcoeff: wavelet coefficients, should be a column vector
%p : resolution of reconstruction (up to 2p)
%j0 : dilation factor of the given wavelet coefficients
function R = get2DHaarReconstruction(wcoeff, p, j0)
	   
    filter = [1 1] ./ sqrt(2);
    sc = IWT2_PO(real(wcoeff),0,filter)+1i* IWT2_PO(imag(wcoeff),0,filter);
    f_recon = @(x,y) sc(floor(x*2^j0)+1,floor(y*2^j0)+1)*2^j0;
    
    R = zeros(2^p);
    x_supp= 0:2^-p:1-2^-p;
    for t = 0:2^p-1
        R(:,t+1) = f_recon(t/2^p,x_supp);
    end
    
end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
