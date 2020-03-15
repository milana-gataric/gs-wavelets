% Set up script for CDJV: Calculate the fourier transforms of the scaling function, right boundary
% scaling functions, and left boundary scaling functions.
%
% Input:
%       a: number of vanishing moments, from the set {2,...,8}
%       eps: Fourier sampling density, require 1/eps to be a natural number
%       S: natural number, this determines the Fourier sampling range: -S*2^R/eps, ... , 2^R*S/eps-1
%       R: maximum scale of wavelet coefficients
% Output:
%       ft_sca: evaluations of the Fourier transform of the main scaling function
%       ft_sca_L: evaluations of the Fourier transform of the left scaling function
%       ft_sca_R: evaluations of the Fourier transform of the right scaling function

function [ft_sca_L,ft_sca, ft_sca_R] = CDJV_Setup(a,eps,S,R)

L= 2^R/eps; 

W = ((S*L:-1:-S*L+1)./L)';

[ft_sca_L, ft_sca, ft_sca_R] = CDJV_FT_Vec(a, W);

end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
