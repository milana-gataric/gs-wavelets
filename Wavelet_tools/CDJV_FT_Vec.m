% Calculate the fourier transforms of the scaling function, right boundary
% scaling functions, and left boundary scaling functions.
%
% Inputs:
%           M :     vector, will calculate the fourier transform at  2*pi*M
%           a :     number of vanishing moments
%
% Ouputs:
%           ft_sca: \int \phi(x) exp(-2\pi M x) dx
%           ft_sca_L(:,j): \int \phi_{L,j}(x) exp(-2\pi M x) dx where \phi_{L,j} is
%           the j^{th} left boundary scaling function
%           ft_sca_R(:,j): \int \phi_{R,j}(x) exp(-2\pi M x) dx where \phi_{R,j} is
%           the j^{th} right boundary scaling function

function [ft_sca_L, ft_sca, ft_sca_R] = CDJV_FT_Vec(a, M)


[LPF,LLPEF,RLPEF] = MakeCDJVFilter('LowPass',a);


ft_sca = zeros(length(M),1);


%%compute the Fourier transforms of the left and right scaling functions
%%evaluated at zero.
ft_0_L = (eye(a)-(1/sqrt(2))*LLPEF(:,1:a))\(LLPEF(:,a+1:end)*(1/sqrt(2))*ones(2*a-1,1));
ft_0_R = (eye(a)-(1/sqrt(2))*RLPEF(:,1:a))\(RLPEF(:,a+1:end)*(1/sqrt(2))*ones(2*a-1,1));


ft_sca_R = zeros(length(M), a);
ft_sca_L = zeros(length(M), a);
for k=1:length(M)
    ft_sca_L(k, :) = CDJV_FTfromFilter(LLPEF, LPF, ft_0_L,2*pi*M(k),15, 0);
    ft_sca_R(k, :) = CDJV_FTfromFilter(RLPEF, LPF, ft_0_R,2*pi*M(k),15, 1);
    ft_sca(k) = FTfromFilter(LPF, 2*pi*M(k),70, 1);
end

ft_sca_R = fliplr(ft_sca_R);


end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
