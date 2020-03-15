%% Computes the GS reconstruction matrix for 2D Haar wavelets and a uniform 
%  sampling scheme

clear all;

%% Input Parameters (change as needed)

eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 4; %maximum scale of wavelet coefficients

L=2^R/eps;

%% Computation of the GS matrix 

W=(eps/2^R)*(-1*S*L:S*L-1)';
ft_sca= haar_phi_ft(-1*W);


filter = [1 1] ./ sqrt(2);


iden_R = eye(2^R);
U=zeros(S*L*2,2^R);

for index = 1:2^R
    wc = iden_R(:,index);
    
    y_concat = Haar_Op_Handle(1, wc, S, eps, R, ft_sca);
    
    U(:,index ) = y_concat;

end

%analysis
smallest_sing_v = sqrt(min(eig(U'*U)));
display(smallest_sing_v)
figure();  subplot(1, 2, 1);
imshow(real(U'*U))
subplot(1, 2, 2);
imshow(sqrt(abs(U)));


%% FOR THE TRANSPOSE
% 
% iden_M = eye(2*S*L);
% U_t=zeros(2^R, S*L*2);
% 
% for index = 1:2*S*L
%     wc = iden_M(:,index);
%     
%     y_concat = Haar_Op_Handle(0, wc, S, eps, R, ft_sca);
%     
%     U_t(:,index ) = y_concat;
% 
% end






% Copyright (c) 2014. Clarice Poon and Milana Gataric
