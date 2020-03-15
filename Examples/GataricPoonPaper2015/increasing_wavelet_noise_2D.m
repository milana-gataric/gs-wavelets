%  2d wavelet reconstruction of a CONTINUOUS function when Fourier samples 
%  are taken uniformly (same example is demonstrated for nonuniform samples)
%  the example produces GS reconstruction for different reconstruction
%  basis (Haar, DB2, DB3) and the truncated Fourier reperentation, 
%  with and without noise

%  comment: show superior GS reconstruction and adventages of basis change
%  even with small number of reconstruction vectors (this makes TF obiously bad)

clear all; 
close all;

p=12;

f = @(x,y) sin(5*pi*x).*cos(3*pi*y);

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);
orig=f(x_supp,y_supp).';

figure('Name','Original');
imshow(real(orig)); 

figure('Position',[50,50,1400,700],'Name','Reconstructions')

%% Input Parameters Uniform sampling

R = 6; %maximum scale of wavelet coefficients

eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1

M = S*2^R/eps;

%%  Fourier samples

disp('Calculating the Fourier samples...')
omega = Samples(eps, M, 'sin(5*pi*x)*cos(3*pi*y)', 0, 1, 0, 1 );

%% Haar wavelets

disp('Haar wavelets')
ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');

%% Compute the Generalized function GS_handle and solve for wavelet coefficients z

disp('Computing Haar reconstruction...');
GS_handle = @(x,mode) Haar_Op2_VecHandle(mode, x, S, eps, R, ft_sca);
z = lsqr(GS_handle,omega(:),10^(-12));
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z))).';

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
reconim = get2DHaarReconstruction(wcoeff, p, R); 

%% Calulate error of GS reconstruction 

error_gs_haar=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_haar); disp(' ');

subplot(2,4,1); imshow(real(reconim)); title({'GS Haar'; ['Error: ',num2str(error_gs_haar)]}); 

%% DB2 wavelets

disp('DB2 wavelets')
a = 2; %number of vanishing moments

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Computing Fourier transform of the scaling functions...')
[ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);

%% Compute the Generalized function GS_handle and solve for wavelet coefficients z

disp('Computing DB2 reconstruction...');
GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
z = lsqr(GS_handle,omega(:),10^(-12));
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z))).';

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
reconim = get2DReconstruction(wcoeff, a, p, R);

%% Calulate error of GS reconstruction 

error_gs_db2=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_db2); disp(' ');

subplot(2,4,2); imshow(real(reconim)); title({'GS DB2'; ['Error: ',num2str(error_gs_db2)]}); 

%% DB3 wavelets

disp('DB3 wavelets')
a = 3; %number of vanishing moments

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Computing Fourier transform of the scaling functions...')
[ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);

%% Compute the Generalized function GS_handle and solve for wavelet coefficients z

disp('Computing DB3 reconstruction...');
GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
z = lsqr(GS_handle,omega(:),10^(-12));
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z))).';

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
reconim = get2DReconstruction(wcoeff, a, p, R);

%% Calulate error of GS reconstruction 

error_gs_db3=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_db3); disp(' ');

subplot(2,4,3); imshow(real(reconim)); title({'GS DB3'; ['Error: ',num2str(error_gs_db3)]}); 

%% Truncated Fourier representation

disp('Computing truncated Fourier representation...')
f_pad = zeros(2*M/eps);
omega_sft = fftshift(omega);
f_pad(1:M,1:M) = omega_sft(1:M,1:M);
f_pad(end - M+1:end, end -M+1:end) = omega_sft(M+1:2*M,M+1:2*M);
f_pad(1:M, end -M+1:end) = omega_sft(1:M,M+1:2*M);
f_pad(end -M+1:end, 1:M) = omega_sft(M+1:2*M, 1:M);
fft_im = eps* fft2(f_pad);
fft_rec = fft_im(1:2*M, 1:2*M);

r=2^p/(2*M);
if r>1
    % producing a higher resolution image
    fft_rec_rez=zeros(2^p,2^p);
    for j=1:2*M
        for i=1:2*M
            fft_rec_rez((j-1)*r+1:j*r,(i-1)*r+1:i*r)=fft_rec(j,i)*ones(r,r);
        end
    end
    fft_rec=fft_rec_rez;
elseif r<1
    fft_rec_rez=zeros(2^p,2^p);
    for j=1:2^p
        for i=1:2^p
            fft_rec_rez(i,j) = fft_rec(1/r*i-1,1/r*j-1);
        end
    end
    fft_rec=fft_rec_rez;
end

error_fft = sqrt(sum(sum(abs(orig - fft_rec).^2))/2^(2*p));
fprintf('Error of the Fourier truncation is %d \n', error_fft); disp(' ');

subplot(2,4,4); imshow(real(fft_rec)); title({'Truncated Fourier series'; ['Error: ',num2str(error_fft)]}); 

%% Add noise and repeat

disp('Reconstructions from noisy measurements...')

SNR=30;
noise = randn(size(omega)) + 1i*randn(size(omega));
noise = 1/SNR*norm(omega(:))/norm(noise(:))*noise;
omega_n = omega + noise;

%% Haar wavelets

disp('Haar wavelets')
ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');

%% Compute the Generalized function GS_handle and solve for wavelet coefficients z

disp('Computing Haar reconstruction...');
GS_handle = @(x,mode) Haar_Op2_VecHandle(mode, x, S, eps, R, ft_sca);
z = lsqr(GS_handle,omega_n(:),10^(-12));
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z))).';

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
reconim = get2DHaarReconstruction(wcoeff, p, R); 

%% Calulate error of GS reconstruction 

error_gs_haar=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_haar); disp(' ');

subplot(2,4,5); imshow(real(reconim)); title({'GS Haar'; ['Error: ',num2str(error_gs_haar)]}); 

%% DB2 wavelets

disp('DB2 wavelets')
a = 2; %number of vanishing moments

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Computing Fourier transform of the scaling functions...')
[ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);

%% Compute the Generalized function GS_handle and solve for wavelet coefficients z

disp('Computing DB2 reconstruction...');
GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
z = lsqr(GS_handle,omega_n(:),10^(-12));
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z))).';

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
reconim = get2DReconstruction(wcoeff, a, p, R);

%% Calulate error of GS reconstruction 

error_gs_db2=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_db2); disp(' ');

subplot(2,4,6); imshow(real(reconim)); title({'GS DB2'; ['Error: ',num2str(error_gs_db2)]}); 

%% DB3 wavelets

disp('DB3 wavelets')
a = 3; %number of vanishing moments

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Computing Fourier transform of the scaling functions...')
[ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);

%% Compute the Generalized function GS_handle and solve for wavelet coefficients z

disp('Computing DB3 reconstruction...');
GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
z = lsqr(GS_handle,omega_n(:),10^(-12));
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z))).';

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
reconim = get2DReconstruction(wcoeff, a, p, R);

%% Calulate error of GS reconstruction 

error_gs_db3=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_db3); disp(' ');

subplot(2,4,7); imshow(real(reconim)); title({'GS DB3'; ['Error: ',num2str(error_gs_db3)]}); 

%% Truncated Fourier representation

disp('Computing truncated Fourier representation...')
f_pad = zeros(2*M/eps);
omega_sft = fftshift(omega_n);
f_pad(1:M,1:M) = omega_sft(1:M,1:M);
f_pad(end - M+1:end, end -M+1:end) = omega_sft(M+1:2*M,M+1:2*M);
f_pad(1:M, end -M+1:end) = omega_sft(1:M,M+1:2*M);
f_pad(end -M+1:end, 1:M) = omega_sft(M+1:2*M, 1:M);
fft_im = eps* fft2(f_pad);
fft_rec = fft_im(1:2*M, 1:2*M);

r=2^p/(2*M);
if r>1
    % producing a higher resolution image
    fft_rec_rez=zeros(2^p,2^p);
    for j=1:2*M
        for i=1:2*M
            fft_rec_rez((j-1)*r+1:j*r,(i-1)*r+1:i*r)=fft_rec(j,i)*ones(r,r);
        end
    end
    fft_rec=fft_rec_rez;
elseif r<1
    fft_rec_rez=zeros(2^p,2^p);
    for j=1:2^p
        for i=1:2^p
            fft_rec_rez(i,j) = fft_rec(1/r*i-1,1/r*j-1);
        end
    end
    fft_rec=fft_rec_rez;
end

error_fft = sqrt(sum(sum(abs(orig - fft_rec).^2))/2^(2*p));
fprintf('Error of the Fourier truncation is %d \n', error_fft)

subplot(2,4,8); imshow(real(fft_rec)); title({'Truncated Fourier series'; ['Error: ',num2str(error_fft)]}); 


% Copyright (c) 2014. Clarice Poon and Milana Gataric
