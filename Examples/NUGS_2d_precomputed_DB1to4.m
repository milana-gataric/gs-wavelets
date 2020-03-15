%% 2d wavelet reconstruction when Fourier samples are taken on a polar sampling scheme; 
%  file produces the GS reconstuction and the gridding reconstruction

%  The file uses precomputed data (available from http://www.damtp.cam.ac.uk/research/afha/code/).

clear all;
close all;

%% Precomputed data is for the following parametars: 
% 
% D = sqrt(2)/4; % the bound for the density in the Euclidean norm i.e. delta < D
% ss = 'polar'; % type of sampling scheme
% 
% a \in {1,2,3,4} % number of vanishing moments
R = 6; % maximum scale of wavelet coefficients
 
N = 2^R; % number of reconstruction vectors
K = N; % band-width of the samples, i.e. samples are taken from the set Y_K a subset of the ball with radius K

f = @(x,y) sin(5*pi*x).*cos(3*pi*y); % function to be reconstruced

%% Set rezolution

p=12; 

%% Load precomputed sampling scheme with density compensation factors:

load('polar_ss_64');

%% Compute samples from the high rezolution image of the function f(x,y):

% [Xs,Ys] = meshgrid(0+0.5/10000:1/10000:1, 0+0.5/10000:1/10000:1);
% ffs = f(Xs,Ys).';
% S=length(ffs);
% sts = nufft_init(2*pi*1/S*sp,[S S],[10 10],[2*S 2*S],[0 0]);
% omega_00 = nufft((1/S)^2*ffs,sts);

load('samples_papercont_pol64');
omega = sqrt(mu).*omega_00;

%% initialize NUFFT 

disp('Data loaded. Initializing NUFFT...')
stx = nufft_init(2*pi*1/2^R*sp(:,1),2^R,10,2*2^R,0,'minmax:tuned');
sty = nufft_init(2*pi*1/2^R*sp(:,2),2^R,10,2*2^R,0,'minmax:tuned');
st = nufft_init(2*pi*1/2^R*sp,[2^R 2^R],[10 10],[2*2^R 2*2^R],[0 0],'minmax:tuned');


%% Haar

disp('Computing Haar reconstruction...');
ft_sca_x = haar_phi_ft(sp(:,1)/2^R); ft_sca_y = haar_phi_ft(sp(:,2)/2^R);
GS_handle_Haar = @(x,mode) Haar_Op2_VecHandle_NU(x, mode, R, ft_sca_x, ft_sca_y, mu, st);
z_Haar = lsqr(GS_handle_Haar,omega(:),10^(-12),50);
wcoeff_Haar = reshape(z_Haar, sqrt(length(z_Haar)), sqrt(length(z_Haar))).';
reconim_Haar = get2DHaarReconstruction(wcoeff_Haar, p, R);

%% DB2

disp('Computing DB2 reconstruction...');
load('ft_sca_64DB2_polar64')
GS_handle_DB2 = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, 2, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
z_DB2 = lsqr(GS_handle_DB2,omega(:),10^(-12),50);
wcoeff_DB2 = reshape(z_DB2, sqrt(length(z_DB2)), sqrt(length(z_DB2))).';
reconim_DB2 = get2DReconstruction(wcoeff_DB2, 2, p, R);

%% DB3

disp('Computing DB3 reconstruction...');
load('ft_sca_64DB3_polar64')
GS_handle_DB3 = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, 3, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
z_DB3 = lsqr(GS_handle_DB3,omega(:),10^(-12),50);
wcoeff_DB3 = reshape(z_DB3, sqrt(length(z_DB3)), sqrt(length(z_DB3))).';
reconim_DB3 = get2DReconstruction(wcoeff_DB3, 3, p, R);

%% DB4

disp('Computing DB4 reconstruction...');
load('ft_sca_64DB4_polar64')
GS_handle_DB4 = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, 4, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
z_DB4 = lsqr(GS_handle_DB4,omega(:),10^(-12),50);
wcoeff_DB4 = reshape(z_DB4, sqrt(length(z_DB4)), sqrt(length(z_DB4))).';
reconim_DB4 = get2DReconstruction(wcoeff_DB4, 4, p, R);

%% Calculate error 

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);    
orig=f(x_supp,y_supp).';

error_Haar=sqrt(sum(sum(abs(real(orig - reconim_Haar)).^2))/2^(2*p));
fprintf('Haar GS reconstruction error is %d \n',error_Haar); 

error_DB2=sqrt(sum(sum(abs(real(orig - reconim_DB2)).^2))/2^(2*p));
fprintf('DB2 GS reconstruction error is %d \n',error_DB2); 

error_DB3=sqrt(sum(sum(abs(real(orig - reconim_DB3)).^2))/2^(2*p));
fprintf('DB3 GS reconstruction error is %d \n',error_DB3); 

error_DB4=sqrt(sum(sum(abs(real(orig - reconim_DB4)).^2))/2^(2*p));
fprintf('DB3 GS reconstruction error is %d \n',error_DB4); 

%% Show images 

figure('Position', [100, 100, 1350, 500],'Name','GS Reconstruction of 2D continious function from polar samples');
subplot(1,4,1); imshow(real(reconim_Haar)); title({'Haar GS reconstruction'; ['Error: ',num2str(error_Haar)]})
subplot(1,4,2); imshow(real(reconim_DB2)); title({'DB2 GS reconstruction'; ['Error: ',num2str(error_DB2)]})
subplot(1,4,3); imshow(real(reconim_DB3)); title({'DB3 GS reconstruction'; ['Error: ',num2str(error_DB3)]})
subplot(1,4,4); imshow(real(reconim_DB4)); title({'DB4 GS reconstruction'; ['Error: ',num2str(error_DB4)]})


% Copyright (c) 2014. Clarice Poon and Milana Gataric
