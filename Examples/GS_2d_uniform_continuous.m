% 2d wavelet reconstruction of a continuous function when Fourier samples are taken uniformly;
% the example produces GS reconstruction and truncated Fourier
% representation, it also plots wavelet coefficients;
% full sampling is used; to use sparsity change the file appropriately (or see GSCS_2D.m)

clear all; 
close all;

%% Input Parameters (change as needed)

a = 4; %number of vanishing moments; it can be in {1,2,3,4,5,6,7,8}
eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 9; %maximum scale of wavelet coefficients

M = S*2^R/eps;

f=@(x,y) (sin(10*x).*cos(30*y)).*x; %function to be reconstructed (if changed, samples need to be changed too!)

fprintf('Using %d^2 DB%d wavelets and %d^2 Fourier samples for GS reconstruction...\n',2^R,a,2*M);

tTotal_start = tic;

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Precomputing Fourier transforms of the scaling functions...')
if a==1  % the case of Haar wavelets
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');
else     % the case of boundary corrected wavelets
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);
end

%% mask for the Fourier coefficients; default: full sampling; 

disp('Using full sampling...')
mask = ones((2*M)); % full sampling

% disp('Taking sparse samples...')
% p1=2;p2=3; 
% [mask,sparsity] = Sampling_map(p1,p2,2*M);
% %or
% mask = zeros(M*2);
% mask(M-2^R/4:M+2^R/4,M-2^R/4:M+2^R/4) = 1;
% %or...

%% Compute Fourier samples

disp('Calculating Fourier samples...')

omega = mask.*Samples(eps, M, 'x*(sin(10*x)*cos(30*y))', 0, 1, 0, 1 );

%% Compute the Generalized function GS_handle and solve for wavelet coefficients z
 
disp('Computing GS reconstruction...')
tGS_start = tic;  
if a==1
    GS_handle = @(x,mode) Haar_Op2_VecHandle(mode, x, S, eps, R, ft_sca, mask);
else
    GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
end
z = lsqr(GS_handle,omega(:),10^(-10));

% for sparse sampling use:
% opts = spgSetParms('verbosity',1);
% z = spg_bpdn(GS_handle,omega(:),0.0005,opts);

wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z)));

%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
p=12;
if a==1
    reconim = get2DHaarReconstruction(wcoeff, p, R);
else
    reconim = get2DReconstruction(wcoeff, a, p, R);
end

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);
orig=f(x_supp,y_supp);
figure('Name', 'Original function'); imagesc(real(orig))
figure('Name', 'GS reconstruction'); imagesc(real(reconim)); 

tGS = toc(tGS_start); 
tTotal = toc(tTotal_start); 
fprintf('Total computing time is %d minutes and %f seconds, out of which time for GS is %d minutes and %f seconds \n',floor(tTotal/60),rem(tTotal,60),floor(tGS/60),rem(tGS,60))

%% Calulate error of GS reconstruction

error_gs=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('GS reconstruction error is %d \n',error_gs); 

title({'GS error', num2str(error_gs)});

%% Plot the wavelet coefficients

J=log2(sqrt(length(z)));
if a==1
    j0=0;
else
    j0=ceil(log2(2*a))+1;
end
ind = [2^(2*j0) 3*2.^(2*(j0:J-1))];
w_vec = zeros(sum(ind),1);

w1 = wcoeff(1:2^j0,1:2^j0);
w_vec(1:ind(1)) = w1(:);
for j=j0:J-1
    w1 = wcoeff(1:2^j, 2^j+1:2^(j+1));
    w2 = wcoeff(2^j+1:2^(j+1), 2^j+1:2^(j+1));
    w3 = wcoeff(2^j+1:2^(j+1), 1:2^j);
    w_vec(ind(j-j0+1)+1:ind(j-j0+1)+ind(j-j0+2)) = [w1(:) ;w2(:) ; w3(:)];
end
figure('Name', 'Reconstructed wavelet coefficients'); stem(abs(w_vec)); set(gca,'xscal','log')

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
    % producing a higher resolution image
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
figure('Name', 'Truncated Fourier series'); imagesc(real(fft_rec));
title({'TF error', num2str(error_fft)});

%%

clear('j','ind','w1','w2','w3','omega_sft','r', 'i')



% Copyright (c) 2014. Clarice Poon and Milana Gataric
