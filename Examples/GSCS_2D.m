%% Compressed sensing with a GS matrix in 2D (uses spgl1)

% The example produces GS reconstruction, truncated Fourier representation 
% of a discontinuous funtion when Fourier samples are taken uniformly,
% and makes 1d intersections for the purpose of comparison.
% Two sampling masks are used: 1. multilevel random sampling map, and 2. sampling map which indexes only samples of low frequencies.

clear all; 
close all;

%% Input Parameters (change as needed)

a = 4; %number of vanishing moments; it can be in {1,2,3,4,5,6,7}
eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 9; %maximum scale of wavelet coefficients
sparsity = 0.05;
M = S*2^R/eps;

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Precomputing Fourier transforms of the scaling functions...')
if a==1  % the case of Haar wavelets
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');
else     % the case of boundary corrected wavelets
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);
end

%% create two masks for the Fourier coefficients
p1=2;p2=4.5; p3=10;
mask = GetSampleMask2D(2*M, ceil(sparsity*M^2*4), p1,p2,p3,100);
fprintf('Taking %.2f%% of the samples...\n',sum(mask(:))*100/M^2/4)

samp_sum = floor(sqrt(sum(mask(:)))/2);
mask_centre =zeros(2*M);
mask_centre(M-samp_sum:M+samp_sum,M-samp_sum:M+samp_sum) = 1;

figure('Position', [150, 150, 900, 400],...
    'Name',sprintf('Sampling masks indexing %d% of the Fourier samples with frequencies up to %d', sparsity*100, S*2^R));
subplot(1,2,1); imagesc(-M:M-1, -M:M-1, mask); title('Mask 1'); 
subplot(1,2,2); imagesc(-M:M-1, -M:M-1, mask_centre); title('Mask 2'); 


%% Define a function and its Fourier samples

disp('Calculating the Fourier samples...')

h = @(x,y) sin(5*pi*x).*cos(3*pi*y) + (x>0.3).*(x<0.58).*(y>0.15).*(y<0.5);
omega_main = Samples(eps, M, 'sin(5*pi*x)*cos(3*pi*y)', 0, 1, 0, 1 ) + Samples(eps, M, '1', 0.3, 0.58, 0.15, 0.5);

slen = 1.5/2^R;
g = @(x,y) (y>0.3).*(y<0.5).*(x>0.5).*(x<0.5+slen)+...
    (y>0.3).*(y<0.5).*(x>0.5+2*slen).*(x<0.5+3*slen)...
    +(y>0.3).*(y<0.5).*(x>0.5+4*slen).*(x<0.5+5*slen)...
    +(y>0.3).*(y<0.5).*(x>0.5+6*slen).*(x<0.5+7*slen);
omega_details = Samples(eps, M, '1', 0.5, 0.5+slen, 0.3, 0.5)+...
    Samples(eps, M, '1', 0.5+2*slen, 0.5+3*slen, 0.3, 0.5)+...
    Samples(eps, M, '1', 0.5+4*slen, 0.5+5*slen, 0.3, 0.5)+...
    Samples(eps, M, '1', 0.5+6*slen, 0.5+7*slen, 0.3, 0.5);

f=@(x,y) h(x,y)+g(x,y);
omega_full = omega_main + omega_details; clearvars omega_main omega_detailsp
%% Compute the Generalized function GS_handle and solve for wavelet coefficients with mask
 
disp('Computing GS reconstruction...')
omega = mask.*omega_full;
vec = @(x) x(:);
if a==1
    GS_handle = @(x,mode) Haar_Op2_VecHandle(mode, x, S, eps, R, ft_sca, mask);
else
    GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
end
opts = spgSetParms('verbosity',1, 'iterations', 50);
z = spg_bpdn(GS_handle,vec(omega),0.0005,opts);
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z)));
clearvars z GS_handle;


%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
p=10;
if a==1   
    reconim = get2DHaarReconstruction(wcoeff, p, R);    
else
    reconim = get2DReconstruction(wcoeff, a, p, R);
end

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);

orig=f(x_supp,y_supp); clearvars x_supp y_supp;
figure('Position', [150, 150, 1250, 500],'Name','2D reconstruction of a discontinuous function using sampling mask 1');
subplot(1,3,1); imshow(real(orig)); title('Original function'); 
subplot(1,3,2); imshow(real(reconim)); title('GS reconstruction'); 


%% Calulate error of GS reconstruction

error_gs=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('GS reconstruction error is %d \n',error_gs);


%% Truncated Fourier representation

disp('Computing truncated Fourier representation...')
omega_sft = fftshift(omega);
fft_im = eps* fft2(omega_sft); clearvars omega_sft;
fft_rec = fft_im(1:2*M, 1:2*M).'; clearvars fft_im;

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
clearvars fft_rec_rez;

error_fft = sqrt(sum(sum(abs(orig - fft_rec).^2))/2^(2*p));
fprintf('Error of the Fourier truncation is %d \n', error_fft)
subplot(1,3,3); imshow(real(fft_rec)); title('Truncated Fourier series'); 

%% Plot 1D intersections
xline=ceil(0.4*2^p);
figure('Position', [150, 150, 1250, 300],'Name','1D intersections of a 2D function (mask 1)');
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

%% Now repeat for the second mask, mask_centre
disp('Computing reconstructions from the second mask...')
omega = mask_centre.*omega_full;

%% Compute the Generalized function GS_handle and solve for wavelet coefficients with mask_centre
disp('Computing GS reconstruction...')

if a==1
    GS_handle = @(x,mode) Haar_Op2_VecHandle(mode, x, S, eps, R, ft_sca, mask_centre);
else
    GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask_centre);
end
opts = spgSetParms('verbosity',1, 'iterations', 50);
z = spg_bpdn(GS_handle,vec(omega),0.0005,opts);
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z))); clearvars GS_handle z
 
%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
p=10;
if a==1    
    reconim = get2DHaarReconstruction(wcoeff, p, R);    
else
    reconim = get2DReconstruction(wcoeff, a, p, R);
end

figure('Position', [150, 150, 1250, 500],'Name','2D reconstruction of a discontinuous function using sampling mask 2');
subplot(1,3,1); imshow(real(orig)); title('Original function'); 
subplot(1,3,2); imshow(real(reconim)); title('GS reconstruction'); 

%% Calulate error of GS reconstruction

error_gs_2=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('GS reconstruction error is %d \n',error_gs_2);

%% Truncated Fourier representation

disp('Computing truncated Fourier representation...')
omega_sft = fftshift(omega);
fft_im = eps* fft2(omega_sft);
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
elseif r<1
    fft_rec_rez=zeros(2^p,2^p);
    for j=1:2^p
        for i=1:2^p
            fft_rec_rez(i,j) = fft_rec(1/r*i-1,1/r*j-1);
        end
    end
    fft_rec=fft_rec_rez;
end
clearvars fft_rec_rez;

error_fft_2 = sqrt(sum(sum(abs(orig - fft_rec).^2))/2^(2*p));
fprintf('Error of the Fourier truncation is %d \n', error_fft_2)
subplot(1,3,3); imshow(real(fft_rec)); title('Truncated Fourier series'); 

%% Plot 1D intersections
xline=ceil(0.4*2^p);
figure('Position', [150, 150, 1250, 300],'Name','1D intersections of a 2D function (mask 2)');
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

clear('j','ind','w1','w2','w3','omega_sft','r', 'i','Ymin', 'Ymax','tGS_start','tTotal_start')

% Copyright (c) 2014. Clarice Poon and Milana Gataric
