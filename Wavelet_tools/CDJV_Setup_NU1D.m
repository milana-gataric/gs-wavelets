% Set up script for CDJV: Calculate the fourier transforms of the scaling function,
% right boundary scaling functions, and left boundary scaling functions.
%
% Input:
%       a:  number of vanishing moments
%       R:  maximum scale of wavelet coefficients
%       sp: a 1D vector of sampling points
% Output:
%       ft_sca:    evaluations of the Fourier transform of the main scaling function
%       ft_sca_L:  evaluations of the Fourier transform of the left scaling function
%       ft_sca_R:  evaluations of the Fourier transform of the right scaling function


function [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup_NU1D(a, R, sp)

W = 1/(2^R)*sp;

[ft_sca_L, ft_sca, ft_sca_R] = CDJV_FT_Vec(a, W);

end


% Copyright (c) 2014. Clarice Poon and Milana Gataric
