%% 2d wavelet reconstruction when Fourier samples are taken on a polar sampling scheme; 
%  file produces the GS reconstuction and the gridding reconstruction

%  The file uses precomputed data (available from http://www.damtp.cam.ac.uk/research/afha/code/).

clear all;
close all;

%% Precomputed data is for the following parametars: 

% D = sqrt(2)/4; % the bound for the density in the Euclidean norm i.e. delta < D
% ss = 'polar'; % type of sampling scheme

% a = 2; % number of vanishing moments
% R = 7; % maximum scale of wavelet coefficients

% N = 2^R; % number of reconstruction vectors
% K = N; % band-width of the samples, i.e. samples are taken from the set Y_K a subset of the ball with radius K

% f = @(x,y) (sin(3*pi*(y+1/2))+(x+1/4).*cos(2*pi*x.*(y+1/2))).^2; % function to be reconstructed

%% Choose whether you want noisy data (SNR=40) or not (SNR=0):

SNR = 00;

%% Load precomputed sampling scheme with density compensation factors:

load('polar_ss_128');

%% Load precomputed samples from the high rezolution image of the function f(x,y):

f = @(x,y) (sin(3*pi*(y+1/2))+(x+1/4).*cos(2*pi*x.*(y+1/2))).^2;
if SNR==0
    load('samples_cont_orig_pol128');
elseif SNR==40
    load('samples_cont_orig_pol128_40SNR');
end
omega = sqrt(mu).*omega_plain;

%% Load precomputed the fourier transforms of the scaling functions:

load('ft_sca_128DB2_polar128')

%% initialize NUFFT 

disp('Data loaded. Initializing NUFFT...')
stx = nufft_init(2*pi*1/2^R*sp(:,1),2^R,7,2*2^R,0,'minmax:tuned');
sty = nufft_init(2*pi*1/2^R*sp(:,2),2^R,7,2*2^R,0,'minmax:tuned');
st = nufft_init(2*pi*1/2^R*sp,[2^R 2^R],[7 7],[2*2^R 2*2^R],[0 0],'minmax:tuned');

%% Compute the Generalized function GS_handle and solve for wavelet coefficient z

disp('Computing GS reconstruction...')
tGS_start = tic;  
GS_handle = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
z = lsqr(GS_handle,omega(:),10^(-12),50);
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z)));

%% reconstruction from wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
p=12;
reconim = get2DReconstruction(wcoeff, a, p, R);

tGS = toc(tGS_start); 
fprintf('Time for computing GS reconstruction is %d minutes and %f seconds \n',floor(tGS/60),rem(tGS,60))

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);    
orig=f(x_supp,y_supp);
figure('Position', [100, 100, 1350, 900],'Name','Reconstruction of 2D continious function from polar samples');
subplot(2,3,1); imshow(real(orig)); title('Original')
subplot(2,3,4); imshow(real(orig(700:2000,1:1300))); title('Zoom');

subplot(2,3,2); imshow(real(reconim)); title('GS reconstruction')
subplot(2,3,5); imshow(real(reconim(700:2000,1:1300))); title('Zoom');

%% Calulate error of GS reconstruction 

dif_gs=abs(real(orig - reconim));
error_gs=sqrt(sum(sum(dif_gs.^2))/2^(2*p));
fprintf('GS reconstruction error is %d \n',error_gs); 

%% Compute the gridding reconstruction (truncated weighted Fourier series)

disp('Computing Gridding reconstruction...')
stg=nufft_init(2*pi*(1/2^p)*sp,[2^p 2^p],[7 7],[2*2^p 2*2^p],[0 0]);

tGR_start = tic; 
f_gr=nufft_adj(1/4*mu.*omega_plain,stg);
tGR = toc(tGR_start); 
fprintf('Time for computing Gridding reconstruction is %d minutes and %f seconds \n',floor(tGR/60),rem(tGR,60))

subplot(2,3,3); imshow(real(f_gr)); title('Gridding reconstruction');
subplot(2,3,6); imshow(real(f_gr(700:2000,1:1300))); title('Zoom');

dif_gr=abs(orig - f_gr);
error_gr=sqrt(sum(sum(dif_gr.^2))/2^(2*p));
fprintf('Error of Gridding reconstruction is %d \n', error_gr)

%% Plot differences

figure('Position', [200, 200, 950, 400],'Name','Difference between reconstruction and original');
if SNR==0
    subplot(1,2,1); imshow(dif_gs,[0,0.003]); title('GS difference from original')
    colorbar
    subplot(1,2,2); imshow(dif_gr,[0,0.05]); title('Gridding difference from original')
    colorbar
elseif SNR==40
    subplot(1,2,1); imshow(dif_gs,[0,0.01]); title('GS difference from original')
    colorbar
    subplot(1,2,2); imshow(dif_gr,[0,0.1]); title('Gridding difference from original')
    colorbar
end


% Copyright (c) 2014. Clarice Poon and Milana Gataric
