%% Generates 1D CDJV matrix when Fourier samples are taken uniformly; 

clear all;

%% Parameters

a = 2; % number of vanishing moments; a can be in {1,2,3}
R = 5; % maximum scale of wavelet coefficients
eps = 1; % distance between consecutive samples; require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; % Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1

L = 2^R/eps;

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

[ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);

%% Computation of the GS matrix

iden_R = eye(2^R);
U=zeros(S*L*2,2^R);

for index = 1:2^R
    wc = iden_R(:,index);

    y_concat = CDJV_Op_Handle(1,wc, S,  eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
    U(:,index ) = y_concat;

end

%analysis
min_sing =sqrt(min(eig(U'*U)));
display(min_sing);

% figure();  subplot( 1, 2, 1);
% imshow(U'*U); title(min_sing);
% subplot(1, 2, 2); 
% imshow(sqrt(abs(U)));

normU = sqrt(max(eig(U'*U)));
display(normU)


%% Generate the transpose matrix directly
%
% iden_M= eye(S*L*2);
% U_t=zeros(2^R, S*L*2);
% for index = 1:S*L*2
%     wc = iden_M(:,index);
% 
%     y_concat = CDJV_Op_Handle(0, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
%     U_t(:,index ) = y_concat;
% 
% end
% 
% figure(); imshow(U_t*U); title(sqrt(min(eig(U_t*U))))


% Copyright (c) 2014. Clarice Poon and Milana Gataric
