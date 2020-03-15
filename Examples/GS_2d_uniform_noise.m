%% 2d wavelet reconstruction when Fourier samples are taken uniformly; 
%  The white Gaussian noise is added to the Fourier samples 

clear all; 
close all;

%% Input Parameters (change as needed)

a=2; %number of vanishing moments; it can be in {1,2,3,4,5,6,7,8}
eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S =1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 7; %maximum scale of wavelet coefficients
SNR = 40; %signal-to-noise ratio

M = S*2^R/eps;

%% Define a function and its Fourier samples

disp('Calculate the Fourier samples and add the noise...')

f=@(x,y) (sin(10*x).*cos(30*y)).*x;
omega_plain = Samples(eps, M, 'x*(sin(10*x)*cos(30*y))', 0, 1, 0, 1 );
noise = randn(size(omega_plain)) + 1i*randn(size(omega_plain));
noise = 1/SNR*norm(omega_plain(:))/norm(noise(:))*noise;
omega = omega_plain + noise;

p=12;
[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);
orig=f(x_supp,y_supp);
figure('Name', 'Original function'); imshow(real(orig))


%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Precomputing Fourier transforms of the scaling functions...')
if a==1  % the case of Haar wavelets
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');
else     % the case of boundary corrected wavelets
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);
end

%% Compute the Generalized function GS_handle and solve for wavelet coefficients z
 
disp('Computing GS reconstruction...')

if a==1
    GS_handle = @(x,mode) Haar_Op2_VecHandle(mode, x, S, eps, R, ft_sca);
else
    GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
end
z = lsqr(GS_handle,omega(:),10^(-12));

wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z)));

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
if a==1
    reconim = get2DHaarReconstruction(wcoeff, p, R);
else
    reconim = get2DReconstruction(wcoeff, a, p, R);
end
figure('Name', 'GS reconstruction'); imshow(real(reconim)); 

%% Calulate error of GS reconstruction 

error_gs=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('GS reconstruction error is %d \n',error_gs); 

title(['GS error: ',num2str(error_gs)]);

%% Truncated Fourier representation

disp('Computing trucated Fourier representation...')
f_pad = zeros(2*M/eps);
omega_sft = fftshift(omega);
f_pad(1:M,1:M) = omega_sft(1:M,1:M);
f_pad(end - M+1:end, end -M+1:end) = omega_sft(M+1:2*M,M+1:2*M);
f_pad(1:M, end -M+1:end) = omega_sft(1:M,M+1:2*M);
f_pad(end -M+1:end, 1:M) = omega_sft(M+1:2*M, 1:M);
fft_im = eps* fft2(f_pad);
fft_rec = fft_im(1:2*M, 1:2*M).';

r=2^p/(2*M);
if r>1
    % producing a higher rezolution image
    fft_rec_rez=zeros(2^p,2^p);

    for j=1:2*M
        for i=1:2*M
            fft_rec_rez((j-1)*r+1:j*r,(i-1)*r+1:i*r)=fft_rec(j,i)*ones(r,r);
        end
    end
    fft_rec=fft_rec_rez;
end

error_fft = sqrt(sum(sum(abs(orig - fft_rec).^2))/2^(2*p));
fprintf('Error of the Fourier truncation is %d \n', error_fft)
figure('Name', 'Truncated Fourier series'); imshow(real(fft_rec));
title(['TF error: ', num2str(error_fft)]);

%%

clear('j','ind','w1','w2','w3','omega_sft','r', 'i','f_pad')

% Copyright (c) 2014. Clarice Poon and Milana Gataric
