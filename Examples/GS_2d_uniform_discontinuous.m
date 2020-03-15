% 2d wavelet reconstruction of a discontinuous function when Fourier samples are taken uniformly; 
% the example produces GS reconstruction, trucated Fourier reperentation and makes 1d intersections for the purpose of comparison
% default: full sampling (to use sparsity see GSCS_2D.m)

clear all; 
close all;

%% Input Parameters (change as needed)

a = 4; %number of vanishing moments; it can be in {1,2,3,4,5,6,7,8}
eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 8; %maximum scale of wavelet coefficients

M = S*2^R/eps;

f = @(x,y) sin(5*pi*x).*cos(3*pi*y) + (x>0.3).*(x<0.58).*(y>0.15).*(y<0.5); %function to be reconstructed (if changed, samples need to be changed too!)

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

%% Fourier samples

disp('Calculating the Fourier samples...')

omega = mask.*Samples(eps, M, 'sin(5*pi*x)*cos(3*pi*y)', 0, 1, 0, 1 ) + Samples(eps, M, '1', 0.3, 0.58, 0.15, 0.5);

%% Compute the Generalized function GS_handle and solve for wavelet coefficients z
 
disp('Computing GS reconstruction...')
tGS_start = tic;  
if a==1
    GS_handle = @(x,mode) Haar_Op2_VecHandle(mode, x, S, eps, R, ft_sca, mask);
else
    GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
end
z = lsqr(GS_handle,omega(:),10^(-12));
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z))).';


%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
p=13;
if a==1
    reconim = get2DHaarReconstruction(wcoeff, p, R); 
else
    reconim = get2DReconstruction(wcoeff, a, p, R);
end

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);
orig=f(x_supp,y_supp).';
figure('Position', [150, 150, 1250, 500],'Name','2D reconstruction of a discontinuous funcion from uniform samples');
subplot(1,3,1); imshow(real(orig)); title('Original function'); 
subplot(1,3,2); imshow(real(reconim)); title('GS reconstruction'); 

tGS = toc(tGS_start); 
tTotal = toc(tTotal_start); 
fprintf('Total computing time is %d minutes and %f seconds, out of which time for GS is %d minutes and %f seconds \n',floor(tTotal/60),rem(tTotal,60),floor(tGS/60),rem(tGS,60))

%% Calulate error of GS reconstruction 

error_gs=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('GS reconstruction error is %d \n',error_gs);

%% Truncated Fourier representation

disp('Computing trucated Fourier representation...')
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
    % producing a higher rezolution image
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
subplot(1,3,3); imshow(real(fft_rec)); title('Truncated Fourier series'); 

%% Plot 1D intersections

xline=(1/2)*2^p;
figure('Position', [150, 150, 1250, 300],'Name','1D intersections of a 2D function');
subplot(1,3,1);
plot(0:2^-p:1-2^-p,real(orig(xline,:)));
Ymax=max(orig(xline,:)); Ymin=min(orig(xline,:)); 
set(gca,'YLim',[Ymin-0.4 Ymax+0.4]); set(gca,'XLim',[-0.03 1.03]);
title('1D Intersection of original function')
subplot(1,3,2);
plot(0:2^-p:1-2^-p,real(reconim(xline,:)));
set(gca,'YLim',[Ymin-0.4 Ymax+0.4]); set(gca,'XLim',[-0.03 1.03]);
title('1D Intersection of GS reconstruction')
subplot(1,3,3);
plot(0:2^-p:1-2^-p,real(fft_rec(xline,:)));
set(gca,'YLim',[Ymin-0.4 Ymax+0.4]); set(gca,'XLim',[-0.03 1.03]);
title('1D Intersection of truncated Fourier series')

%%

clear('j','ind','w1','w2','w3','omega_sft','r', 'i','Ymin', 'Ymax','tGS_start','tTotal_start','xline')

% Copyright (c) 2014. Clarice Poon and Milana Gataric
