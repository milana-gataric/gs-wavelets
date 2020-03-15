% 1d wavelet reconstruction when Fourier samples are taken uniformly; 
% produces GS and truncated Fourier representation as well.

clear all;
close all;

%% Parameters

a = 4; % number of vanishing moments; a can be in {1,2,3,4,5,6,7,8}
R = 8; % maximum scale of wavelet coefficients
eps = 1; % distance between consecutive samples; require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; % Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1

M = S*2^R/eps; % the number of sampling points equals to 2*M

f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6); % the function to be reconstructed (if changed, samples need to be changed too)
% %or 
% f = @(x)  (0.5.*(x>=1/3).*(x<=2/3)) + 0.5.*(x>=2/5).*(x<=(2/5+1/300))+(x>=3/5).*(x<=(3/5+1/300));
% % or define your example...

fprintf('Using %d DB%d wavelets and %d Fourier samples for GS reconstruction...\n',2^R,a,2*M);

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

if a==1  % case of Haar wavelets
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');
else     % case of boundary corrected wavelets
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);
end

%% mask for the Fourier coefficients; default is full sampling 

mask = ones((2*M),1); 

%% Compute Fourier samples

omega = mask.*Samples(eps, M, '(x+1)*(x-0.5)*(x-0.25)*(x+3)*(x-0.6)', 0, 1);

% %or
% omega =  mask.*Samples( eps, M, '0.5', 1/3, 2/3 ) + Samples( eps, M, '0.5', 2/5, 2/5+1/300 ) + Samples( eps, M, '1', 3/5, 3/5+1/300 );

%% Compute the Generalized function GS_handle and solve for wavelet coefficient z

if a==1
    GS_handle = @(x,mode) Haar_Op_Handle(mode, x, S, eps, R, ft_sca, mask);
else
    GS_handle = @(x,mode) CDJV_Op_Handle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
end
z = lsqr(GS_handle,omega(1:2*M),10^(-10));

%% reconstruction from wavelet coefficients

p=16;
if a==1
    recons = get1DHaarReconstruction(z, p, R);    
else
    recons = get1DReconstruction(z, a, p, R);
end

orig=f((2^-p:2^-p:1)');
error_gs = sqrt(sum(abs(orig-recons).^2)/2^(p));
fprintf('GS error is %d \n',error_gs)

figure('Position', [300, 400, 1200, 270],'Name','1D reconstruction from uniform Fourier samples');
subplot(1,3,1); plot(2^-p:2^-p:1, orig); title('Original function');axis([-0.02, 1.02 min(orig)-0.2 max(orig)+0.2])
subplot(1,3,2); plot(2^-p:2^-p:1, real(recons)); title({'GS reconstruction'; ['Error: ',num2str(error_gs)]});axis([-0.02, 1.02 min(orig)-0.2 max(orig)+0.2])

%% Truncated Fourier representation

Supp=2^-p:2^-p:1;
y=zeros(size(Supp));
for j=M:-1:-M+1
    y=y+sqrt(eps)*exp(2*eps*1i*j*Supp*pi)*omega(M-j+1);
end

error_ft = sqrt(sum(abs(f((2^-p:2^-p:1))-y).^2)/2^(p));
fprintf('Truncated Fourier error is %d \n',error_ft)

subplot(1,3,3); plot(Supp, real(y)); title({'Truncated Fourier reconstruction'; ['Error: ',num2str(error_ft)]});axis([-0.02, 1.02 min(orig)-0.2 max(orig)+0.2])    



% Copyright (c) 2014. Clarice Poon and Milana Gataric
