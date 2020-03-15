% Set up script for CDJV: Calculate the fourier transforms of the scaling function,
% right boundary scaling functions, and left boundary scaling functions.
%
% Input:
%       a:  number of vanishing moments
%       R:  maximum scale of wavelet coefficients
%       sp: a 2D vector of sampling points (for x and y coordinate)
% Output:
%       ft_sca_x:    evaluations of the Fourier transform of the main
%                    scaling function at sp(:,1) 
%       ft_sca_L_x:  evaluations of the Fourier transform of the left
%                    scaling function at sp(:,1)
%       ft_sca_R_x:  evaluations of the Fourier transform of the right
%                    scaling function at sp(:,1)
%       ft_sca_y:    evaluations of the Fourier transform of the main
%                    scaling function at sp(:,2)
%       ft_sca_L_y:  evaluations of the Fourier transform of the left
%                    scaling function at sp(:,2)
%       ft_sca_R_y:  evaluations of the Fourier transform of the right
%                    scaling function at sp(:,2)


function [ft_sca_L_x, ft_sca_x, ft_sca_R_x, ft_sca_L_y, ft_sca_y, ft_sca_R_y] = CDJV_Setup_NU2D(a, R, sp)

W = 1/(2^R)*sp;

[ft_sca_L_x, ft_sca_x, ft_sca_R_x] = CDJV_FT_Vec(a, W(:,1));
[ft_sca_L_y, ft_sca_y, ft_sca_R_y] = CDJV_FT_Vec(a, W(:,2));

end


% Copyright (c) 2014. Clarice Poon and Milana Gataric
