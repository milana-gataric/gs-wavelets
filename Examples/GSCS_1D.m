%% Compressed sensing with a GS matrix in 1D (uses spgl1)

%  The example produces GS reconstruction and truncated Fourier
%  representation of a discontinuous funtion when Fourier samples are taken uniformly;
%  Two sampling masks are used: 1. multilevel random sampling map, 
%  and 2. sampling map which indexes only samples of low frequencies.

clear all; 
close all;

%% Input Parameters (change as needed)

a = 4; %number of vanishing moments; it can be in {1,2,3,4,5,6,7}
eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 11; %maximum scale of wavelet coefficients

M = S*2^R/eps;

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Precomputing Fourier transforms of the scaling functions...')
if a==1  % the case of Haar wavelets
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');
else     % the case of boundary corrected wavelets
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);
end

%% create two masks for the Fourier coefficients
p1=2;p2=3.4;
M0=M*2;
mask = GetSampleMask1D(M0,ceil(0.1*M0), p1,p2, 30, 200);

sparsity = sum(mask)/(M0);
fprintf('Taking %.2f%% of the samples...\n',sparsity*100)
figure('Position', [350, 350, 500, 450],'Name','Mask 1'); 
plot(mask, '.'); 

samp_sum = floor(sum(mask(:))/2);
mask_centre =zeros(2*M,1);
mask_centre(M-samp_sum:M+samp_sum) = 1;
figure('Position', [350, 350, 500, 450],'Name','Mask 2');
plot(mask_centre, '.');
%% Define a function and its Fourier samples

disp('Calculating the Fourier samples...')

f = @(x) sin(5*pi*x) + (x>0.3).*(x<0.58);
xylims=[-0.02 1.02 -1.5 3.6];
omega_main = Samples(eps, M, 'sin(5*pi*x)', 0, 1) + Samples(eps, M, '1', 0.3, 0.58);

slen = 4/2^R;
g = @(x) (x>0.5).*(x<0.5+slen)+...
    1.5*(x>0.5+2*slen).*(x<0.5+3*slen)...
    +0.5*(x>0.5+4*slen).*(x<0.5+5*slen)...
    +(x>0.5+6*slen).*(x<0.5+7*slen);
omega_details = Samples(eps, M, '1', 0.5, 0.5+slen)+1.5*...
    Samples(eps, M, '1', 0.5+2*slen, 0.5+3*slen) + 0.5*...
    Samples(eps, M, '1', 0.5+4*slen, 0.5+5*slen)+...
    Samples(eps, M, '1', 0.5+6*slen, 0.5+7*slen);

f=@(x) f(x)+g(x);
omega = mask.*(omega_main + omega_details);
%% Compute the Generalized function GS_handle and solve for wavelet coefficients with mask
 
disp('Computing GS reconstruction...')

if a==1
    GS_handle = @(x,mode) Haar_Op_Handle(mode, x, S, eps, R, ft_sca, mask);
else
    GS_handle = @(x,mode) CDJV_Op_Handle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
end
opts = spgSetParms('verbosity',1, 'iterations', 500);
wcoeff = spg_bpdn(GS_handle,omega(:),0.0005,opts);


%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
p=max(12, R+1);
if a==1
    
    reconim = get1DHaarReconstruction(wcoeff, p, R);
    
else
    reconim = get1DReconstruction(wcoeff, a, p, R);
end

x_supp=(0:2^-p:1-2^-p)';

orig=f(x_supp);
figure('Position', [150, 150, 1250, 500],'Name','1D reconstruction of a discontinuous funcion from sparse uniform samples for Mask 1');
subplot(2,3,1); plot(x_supp, real(orig)); axis(xylims); title('Original function'); 
subplot(2,3,4); plot(x_supp, real(orig)); xlim([0.4,0.6]); ylim([1.5,3.7]);title('Zoom'); 

subplot(2,3,2); plot(x_supp, real(reconim)); axis(xylims); title('GS reconstruction'); 
subplot(2,3,5); plot(x_supp, real(reconim)); xlim([0.4,0.6]); ylim([1.5,3.7]);title('Zoom'); 




%% Calulate error of GS reconstruction

error_gs=sqrt(sum(sum(abs(orig - reconim).^2))/2^(p));
fprintf('GS reconstruction error is %d \n',error_gs);


%% Truncated Fourier representation

disp('Computing trucated Fourier representation...')

Supp=2^-p:2^-p:1;
fft_rec=zeros(size(Supp));
for j=M:-1:-M+1
    fft_rec=fft_rec+sqrt(eps)*exp(2*eps*1i*j*Supp*pi)*omega(M-j+1);
end

error_fft = sqrt(sum(abs(f((2^-p:2^-p:1))-fft_rec).^2)/2^(p));
fprintf('Truncated Fourier error is %d \n',error_fft)

subplot(2,3,3); plot(x_supp, real(fft_rec)); axis(xylims); title('Truncated Fourier series'); 
subplot(2,3,6); plot(x_supp, real(fft_rec)); xlim([0.4,0.6]); ylim([1.5,3.7]);title('Zoom'); 


%% Now repeat for the second mask, mask_centre
disp('Computing reconstructions from the second mask...')

%% Compute the Generalized function GS_handle and solve for wavelet coefficients with mask_centre
omega = mask_centre.*(omega_details + omega_main);
disp('Computing GS reconstruction...')

if a==1
    GS_handle = @(x,mode) Haar_Op_Handle(mode, x, S, eps, R, ft_sca, mask_centre);
else
    GS_handle = @(x,mode) CDJV_Op_Handle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask_centre);
end
opts = spgSetParms('verbosity',1, 'iterations', 500);
wcoeff2 = spg_bpdn(GS_handle,omega,0.0005,opts);
 

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
p=max(14,R+2);
if a==1
    
    recon = get1DHaarReconstruction(wcoeff2, p, R);
    
else
    recon = get1DReconstruction(wcoeff2, a, p, R);
end

x_supp =(0:2^-p:1-2^-p)';

orig=f(x_supp);
figure('Position', [150, 150, 1250, 500],'Name','2D reconstruction of a discontinuous funcion from sparse uniform samples for Mask 2');
subplot(2,3,1); plot(x_supp, real(orig)); axis(xylims); title('Original function'); 
subplot(2,3,4); plot(x_supp, real(orig)); xlim([0.4,0.6]); ylim([1.5,3.7]);title('Zoom'); 

subplot(2,3,2); plot(x_supp, real(recon)); axis(xylims); title('GS reconstruction'); 
subplot(2,3,5); plot(x_supp, real(recon)); xlim([0.4,0.6]); ylim([1.5,3.7]);title('Zoom'); 

%% Calulate error of GS reconstruction

error_gs_2=sqrt(sum(sum(abs(orig - recon).^2))/2^(p));
fprintf('GS reconstruction error is %d \n',error_gs_2);

%% Truncated Fourier representation

disp('Computing trucated Fourier representation...')

Supp=2^-p:2^-p:1;
fft_rec2=zeros(size(Supp));
for j=M:-1:-M+1
    fft_rec2=fft_rec2+sqrt(eps)*exp(2*eps*1i*j*Supp*pi)*omega(M-j+1);
end

error_fft_2 = sqrt(sum(abs(f((2^-p:2^-p:1))-fft_rec2).^2)/2^(p));
fprintf('Error of the Fourier truncation is %d \n', error_fft_2)
subplot(2,3,3); plot(Supp, real(fft_rec2)); axis(xylims); title('Truncated Fourier series'); 
subplot(2,3,6); plot(Supp, real(fft_rec2)); xlim([0.4,0.6]); ylim([1.5,3.7]);title('Zoom'); 

%%


% Copyright (c) 2014. Clarice Poon and Milana Gataric
