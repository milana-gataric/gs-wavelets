%% 2d wavelet reconstruction of a smooth function from sampling 5% of its Fourier coefficients (up to frequency M); 
%  the example produces a reconstruction using a GS matrix, and a
%  reconstruction via the more standard compressed sensing matrix DFT DWT^{-1}


clear all; 
close all;

%% Input Parameters (change as needed)

a = 4; %number of vanishing moments; it can be in {1,2,3,4,5,6,7}
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 8; %maximum scale of wavelet coefficients


M = S*2^R;
eps = 1; %Nyquist rate, set to be one for this experiment

%% create mask for the Fourier coefficients (comment/uncomment to try a different mask)

%%Mask 1: full sampling
% mask = ones(M*2);

%%Mask 2: multilevel sampling
%subsampling parmeters
% p1=2; p2=4; p3 = 10;
% sparsity = 0.05;
% mask = GetSampleMask2D(2*M, ceil(sparsity*4*M^2), p1, p2, p3, 100);

%%Mask 3: "radial lines sampling" (samples are on a grid)
mask = GetLinesMask(2*M, 22);

%display the sampling mask
figure('Position', [150, 150, 1250, 500],'Name','Sampling mask');
imshow(mask); title(sprintf('Sampling mask (%.2f%%)', 100*sum(mask(:))/numel(mask)));


%% Define a function and its Fourier samples

disp('Calculating the Fourier samples...')

f = @(x,y) sin(5*pi*y).*cos(3*pi*x).*exp(-x).*exp(-y);
omega = mask.* Samples(eps, M, 'sin(5*pi*y)*cos(3*pi*x)*exp(-x)*exp(-y)', 0, 1, 0, 1 ) ;

%% Precompute the Fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Precomputing Fourier transforms of the scaling functions...')
if a==1  % the case of Haar wavelets
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');
else     % the case of boundary corrected wavelets
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);
end

%% solve for wavelet coefficients with mask using GS/Standard CS
 
disp('Computing GS reconstruction...')

if a==1
    GS_handle = @(x,mode) Haar_Op2_VecHandle(mode, x, S, eps, R, ft_sca, mask);
else
    GS_handle = @(x,mode) CDJV_Op2_VecHandle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
end
opts = spgSetParms('verbosity',1, 'iterations', 100);
z = spg_bpdn(GS_handle,omega(:),0.0005,opts);
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z)));
clearvars z;

disp('Computing standard CS reconstruction...')

vec = @(x) x(:);
opts = spgSetParms('verbosity',1, 'iterations', 100);
discreteCS_handle = @(x,mode) DFT_DWT_Op(mode, x, a, fftshift(mask.'));
z_discrete = spg_bpdn(discreteCS_handle,vec(fftshift(omega.')),0.0005,opts);
wcoeff_discrete = reshape(z_discrete, sqrt(length(z_discrete)), sqrt(length(z_discrete)));
wcoeff_discrete = wcoeff_discrete.';


%% Compute image from reconstructed wavelet coefficients

disp('Computing image from GS reconstructed wavelet coefficients...')
p=11;
if a==1    
    reconim = get2DHaarReconstruction(wcoeff, p, R);    
else
    reconim = get2DReconstruction(wcoeff, a, p, R);
end

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);

orig=f(x_supp,y_supp); clearvars x_supp y_supp;
figure('Position', [150, 150, 1250, 500],'Name','Comparison between compressed sensing using GS and using DFT DWT^{-1}');
subplot(2,3,1); imshow(real(orig)); title('Original function'); 
subplot(2,3,4); imshow(real(orig(1:300,1:300))); title('Zoom (top right)');


subplot(2,3,2); imshow(real(reconim)); title('GS reconstruction'); 
subplot(2,3,5); imshow(real(reconim(1:300,1:300))); title('Zoom (top right)');


disp('Computing image from DFT DWT^{-1} reconstructed wavelet coefficients...')
if a==1    
    reconim_ds = get2DHaarReconstruction(wcoeff_discrete, p, log2(length(wcoeff_discrete)));    
else
    reconim_ds = get2DReconstruction(wcoeff_discrete, a, p, log2(length(wcoeff_discrete)));
end

subplot(2,3,3); imshow(real(reconim_ds)); title('Standard CS reconstruction'); 
subplot(2,3,6); imshow(real(reconim_ds(1:300,1:300))); title('Zoom (top right)');


%% Calulate error of GS and standard CS reconstruction

error_gs=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('GS reconstruction error is %d \n',error_gs);


error_gs=sqrt(sum(sum(abs(orig - reconim_ds).^2))/2^(2*p));
fprintf('Standard CS reconstruction error is %d \n',error_gs);




clear('j','ind','w1','w2','w3','omega_sft','r', 'i','Ymin', 'Ymax','tGS_start','tTotal_start')







% Copyright (c) 2014. Clarice Poon and Milana Gataric