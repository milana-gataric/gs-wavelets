%  2d wavelet reconstruction of a CONTINUOUS function when Fourier samples 
%  are taken on a polar sampling scheme 
%  the example produces GS reconstruction for different reconstruction
%  basis (Haar, DB2, DB3) and the truncated Fourier reperentation (aka gridding),
%  with and without noise

% Uses precomputed data (available from http://www.damtp.cam.ac.uk/research/afha/code/)

clear all; 
close all;

p=12;

f = @(x,y) sin(5*pi*x).*cos(3*pi*y);

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);
orig=f(x_supp,y_supp).';

figure('Name','Original');
imshow(real(orig)); 

figure('Position',[50,50,1400,700],'Name','Reconstructions from nonuniform samples')

%% Input Parameters Uniform sampling

R = 6; %maximum scale of wavelet coefficients
N = 2^R;
K = N;

load('polar_ss_64');

%% initialize NUFFT 

disp('Initializing NUFFT...')
stx = nufft_init(2*pi*1/2^R*sp(:,1),2^R,10,2*2^R,0,'minmax:tuned');
sty = nufft_init(2*pi*1/2^R*sp(:,2),2^R,10,2*2^R,0,'minmax:tuned');
st = nufft_init(2*pi*1/2^R*sp,[2^R 2^R],[10 10],[2*2^R 2*2^R],[0 0],'minmax:tuned');

%%  Fourier samples

disp('Calculating the Fourier samples...')

% [Xs,Ys] = meshgrid(0+0.5/20000:1/20000:1, 0+0.5/20000:1/20000:1);
% ffs = f(Xs,Ys);
% S=length(ffs);
% sts = nufft_init(2*pi*1/S*sp,[S S],[10 10],[2*S 2*S],[0 0]);
% omega_00 = nufft((1/S)^2*ffs,sts);

load('samples_papercont_pol64');

omega = sqrt(mu).*omega_00;

disp(' ');

%% Haar wavelets
%  Compute the Generalized function GS_handle and solve for wavelet coefficients z

disp('Computing Haar reconstruction...');
ft_sca_x = haar_phi_ft(sp(:,1)/2^R); ft_sca_y = haar_phi_ft(sp(:,2)/2^R);
GS_handle_Haar = @(x,mode) Haar_Op2_VecHandle_NU(x, mode, R, ft_sca_x, ft_sca_y, mu, st);
z_Haar = lsqr(GS_handle_Haar,omega(:),10^(-12),50);
wcoeff_Haar = reshape(z_Haar, sqrt(length(z_Haar)), sqrt(length(z_Haar))).';

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
reconim_Haar = get2DHaarReconstruction(wcoeff_Haar, p, R);

%% Calulate error of GS reconstruction 

error_gs_haar=sqrt(sum(sum(abs(orig - reconim_Haar).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_haar); disp(' ');

subplot(2,4,1); imshow(real(reconim_Haar)); title({'GS Haar'; ['Error: ',num2str(error_gs_haar)]}); 

%% DB2 wavelets

% Compute the Generalized function GS_handle and solve for wavelet coefficients z
disp('Computing DB2 reconstruction...');
%[ft_sca_L_x, ft_sca_x, ft_sca_R_x, ft_sca_L_y, ft_sca_y, ft_sca_R_y] = CDJV_Setup_NU2D(2, R, sp);
load('ft_sca_64DB2_polar64')
GS_handle_DB2 = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, 2, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
z_DB2 = lsqr(GS_handle_DB2,omega(:),10^(-12),50);
wcoeff_DB2 = reshape(z_DB2, sqrt(length(z_DB2)), sqrt(length(z_DB2))).';

% Compute image from reconstructed wavelet coefficients
disp('Computing image from reconstructed wavelet coefficients...')
reconim_DB2 = get2DReconstruction(wcoeff_DB2, 2, p, R);

% Calulate error of GS reconstruction 
error_gs_db2=sqrt(sum(sum(abs(orig - reconim_DB2).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_db2); disp(' ');

subplot(2,4,2); imshow(real(reconim_DB2)); title({'GS DB2'; ['Error: ',num2str(error_gs_db2)]}); 

%% DB3 wavelets

% Compute the Generalized function GS_handle and solve for wavelet coefficients z
disp('Computing DB3 reconstruction...');
load('ft_sca_64DB3_polar64')
GS_handle_DB3 = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, 3, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
z_DB3 = lsqr(GS_handle_DB3,omega(:),10^(-12),50);
wcoeff_DB3 = reshape(z_DB3, sqrt(length(z_DB3)), sqrt(length(z_DB3))).';

% Compute image from reconstructed wavelet coefficients
disp('Computing image from reconstructed wavelet coefficients...')
reconim_DB3 = get2DReconstruction(wcoeff_DB3, 3, p, R);

% Calulate error of GS reconstruction 
error_gs_db3=sqrt(sum(sum(abs(orig - reconim_DB3).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_db3); disp(' ');

subplot(2,4,3); imshow(real(reconim_DB3)); title({'GS DB3'; ['Error: ',num2str(error_gs_db3)]}); 

%% Truncated Fourier representation

disp('Computing Gridding reconstruction...')
stg=nufft_init(2*pi*(1/2^p)*sp,[2^p 2^p],[10 10],[2*2^p 2*2^p],[0 0]);
f_gr=nufft_adj(mu.*omega_00,stg).';

dif_gr=abs(orig - f_gr);
error_gr=sqrt(sum(sum(dif_gr.^2))/2^(2*p));
fprintf('Error of Gridding reconstruction is %d \n', error_gr); disp(' ');

subplot(2,4,4); imshow(real(f_gr)); title({'Truncated Fourier series'; ['Error: ',num2str(error_gr)]}); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADD NOISE AND REPEAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Reconstructions from noisy measurements...'); disp(' ');

SNR=30;
noise = randn(size(omega_00)) + 1i*randn(size(omega_00));
noise = 1/SNR*norm(omega_00(:))/norm(noise(:))*noise;
omega_n = (omega_00 + noise);
omega_nm = sqrt(mu).*omega_n;

%% Haar wavelets

% Compute the Generalized function GS_handle and solve for wavelet coefficients z
disp('Computing Haar reconstruction...');
ft_sca_x = haar_phi_ft(sp(:,1)/2^R); ft_sca_y = haar_phi_ft(sp(:,2)/2^R);
GS_handle_Haar = @(x,mode) Haar_Op2_VecHandle_NU(x, mode, R, ft_sca_x, ft_sca_y, mu, st);
z_Haar = lsqr(GS_handle_Haar,omega_nm(:),10^(-12),50);
wcoeff_Haar = reshape(z_Haar, sqrt(length(z_Haar)), sqrt(length(z_Haar))).';

% Compute image from reconstructed wavelet coefficients
disp('Computing image from reconstructed wavelet coefficients...')
reconim_Haar = get2DHaarReconstruction(wcoeff_Haar, p, R);

% Calulate error of GS reconstruction 
error_gs_haar=sqrt(sum(sum(abs(orig - reconim_Haar).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_haar); disp(' ');

subplot(2,4,5); imshow(real(reconim_Haar)); title({'GS Haar'; ['Error: ',num2str(error_gs_haar)]}); 

%% DB2 wavelets

% Compute the Generalized function GS_handle and solve for wavelet coefficients z
disp('Computing DB2 reconstruction...');
load('ft_sca_64DB2_polar64')
GS_handle_DB2 = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, 2, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
z_DB2 = lsqr(GS_handle_DB2,omega_nm(:),10^(-12),50);
wcoeff_DB2 = reshape(z_DB2, sqrt(length(z_DB2)), sqrt(length(z_DB2))).';

% Compute image from reconstructed wavelet coefficients
disp('Computing image from reconstructed wavelet coefficients...')
reconim_DB2 = get2DReconstruction(wcoeff_DB2, 2, p, R);

% Calulate error of GS reconstruction 
error_gs_db2=sqrt(sum(sum(abs(orig - reconim_DB2).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_db2); disp(' ');

subplot(2,4,6); imshow(real(reconim_DB2)); title({'GS DB2'; ['Error: ',num2str(error_gs_db2)]}); 

%% DB3 wavelets

% Compute the Generalized function GS_handle and solve for wavelet coefficients z
disp('Computing DB3 reconstruction...');
load('ft_sca_64DB3_polar64')
GS_handle_DB3 = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, 3, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
z_DB3 = lsqr(GS_handle_DB3,omega_nm(:),10^(-12),50);
wcoeff_DB3 = reshape(z_DB3, sqrt(length(z_DB3)), sqrt(length(z_DB3))).';

% Compute image from reconstructed wavelet coefficients
disp('Computing image from reconstructed wavelet coefficients...')
reconim_DB3 = get2DReconstruction(wcoeff_DB3, 3, p, R);

% Calulate error of GS reconstruction 
error_gs_db3=sqrt(sum(sum(abs(orig - reconim_DB3).^2))/2^(2*p));
fprintf('Error is %d \n',error_gs_db3); disp(' ');

subplot(2,4,7); imshow(real(reconim_DB3)); title({'GS DB3'; ['Error: ',num2str(error_gs_db3)]}); 

%% Truncated Fourier representation

disp('Computing Gridding reconstruction...')
stg=nufft_init(2*pi*(1/2^p)*sp,[2^p 2^p],[7 7],[2*2^p 2*2^p],[0 0]);
f_gr=nufft_adj(mu.*omega_n,stg).';

dif_gr=abs(orig - f_gr);
error_gr=sqrt(sum(sum(dif_gr.^2))/2^(2*p));
fprintf('Error of Gridding reconstruction is %d \n', error_gr); disp(' ');

subplot(2,4,8); imshow(real(f_gr)); title({'Truncated Fourier series'; ['Error: ',num2str(error_gr)]}); 


% Copyright (c) 2014. Clarice Poon and Milana Gataric
