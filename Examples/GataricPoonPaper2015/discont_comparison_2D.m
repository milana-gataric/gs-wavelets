%  2d wavelet reconstruction of a DISCONTINUOUS function when Fourier samples are taken uniformly; 
%  the example produces GS reconstruction, truncated Fourier representation
%  and one dimensional intersections;

clear all; 
close all;

%% Input Parameters 

a = 4; %number of vanishing moments
eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 8; %maximum scale of wavelet coefficients

M = S*2^R/eps;

p=12;

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Precomputing Fourier transforms of the scaling functions...')
if a==1  % the case of Haar wavelets
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');
else     % the case of boundary corrected wavelets
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);
end

%% Define a function and its Fourier samples

disp('Calculating the Fourier samples...')

h = @(x,y) sin(5*pi*y).*cos(3*pi*x) + (y>0.3).*(y<0.58).*(x>0.15).*(x<0.5);
slen = 2/2^R;
g = @(x,y) (y>0.3).*(y<0.58).*(x>0.4).*(x<0.4+slen)+...
    +(y>0.3).*(y<0.58).*(x>0.4+4*slen).*(x<0.4+5*slen);
f=@(x,y) h(x,y)+g(x,y);

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);
orig=f(x_supp,y_supp);
figure('Position', [150, 150, 1250, 700],'Name','2D reconstruction of a discontinuous funcion from uniform samples');
subplot(2,3,1); imshow(real(orig)); title('Original function'); 
subplot(2,3,4); imshow(real(orig)); title('Original function zoom'); zoom(4) 

omega_main = Samples(eps, M, 'sin(5*pi*y)*cos(3*pi*x)', 0, 1, 0, 1 ) + Samples(eps, M, '1',0.15, 0.5, 0.3, 0.58);
omega_details = Samples(eps, M, '1', 0.4, 0.4+slen, 0.3, 0.58)+...
    Samples(eps, M, '1',0.4+4*slen, 0.4+5*slen,0.3, 0.58);

omega = omega_main + omega_details;

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
    reconim = get2DReconstruction(wcoeff, a, p, R, 24);
end

%% Calulate error of GS reconstruction 

error_gs=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('GS reconstruction error is %d \n',error_gs);

subplot(2,3,2); imshow(real(reconim)); title('GS reconstruction'); 
subplot(2,3,5); imshow(real(reconim)); title('GS reconstruction zoom'); zoom(4);

%% Truncated Fourier representation

disp('Computing truncated Fourier representation...')
omega_sft = fftshift(omega);
fft_im = eps* fft2(omega_sft);
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

subplot(2,3,3); imshow(real(fft_rec)); title('Truncated Fourier series'); 
subplot(2,3,6); imshow(real(fft_rec)); title('Truncated Fourier series zoom'); zoom(4);

%% Plot 1D intersections

xline=(1/2-1/16)*2^p;
figure('Position', [150, 150, 1250, 300],'Name','1D intersections of a 2D function from uniform samples');
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



% Copyright (c) 2014. Clarice Poon and Milana Gataric
