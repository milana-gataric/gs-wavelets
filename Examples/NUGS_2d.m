%% 2d wavelet reconstruction when Fourier samples are taken on a nonuniform 
%  sampling scheme: choose radial (polar) or spiral;
%  file produces both GS and gridding reconstruction as well as the plot of wavelet coefficients
%  default: full sampling & lsqr

% Warning: Precomputation steps take longer time for K>=64 (for K=64, 
% about 10 minutes); use precomputed data for bigger examples

clear all;
close all;

%% Parameters

D = sqrt(2)/4; % the bound for the density in the Euclidean norm 
               % i.e. delta < D; max D is sqrt(2)/4
ss = 'polar'; % type of the sampling scheme, it can be a polar or an (interleaving) 
              % spiral: 'polar', 'spiral 1', 'spiral 2', 'spiral 3', ...
               
a = 4; % number of vanishing moments
R = 5; % maximum scale of wavelet coefficients

N = 2^R; % number of reconstruction vectors
K = N; % band-width of the samples

f = @(x,y) ((1-x).*y+sin(y)).^2; % function to be reconstructed (if changed, samples have to be changed as well!)
% %or a discontinuous example:
% f = @(x,y) (x>0.25).*(x<0.75).*(y>0.25).*(y<0.75);

fprintf('Using %d^2 DB%d wavelets and %s Fourier samples with sampling bandwidth %d for GS reconstruction...\n',2^R,a,ss,K);

tTotal_start = tic;

%% Sampling scheme and density compensation factors (Voronoi regions)
if exist('mri_density_comp', 'file') % the same function has different names across the different versions of the package IR
    density_fn = @(x) mri_density_comp(x, 'voronoi');
else
    density_fn = @(x) ir_mri_density_comp(x, 'voronoi');
end
    
disp('Calculating the sampling scheme and density compensation factors...')
if strcmp(ss,'polar')
    sp = polar_sampling(K, D);
    mu = 4*density_fn(sp); % takes time for large reconstructions
elseif strfind(ss, 'spiral')
    [~, leaves] = strtok(ss);
    sp = spiral_sampling(K, D, str2double(leaves));
    mu = 4*density_fn(sp); % takes time for large reconstructions
    mu(mu>5)=5; % a correction is needed since the function 'mri_density_comp' doesn't work perfectly well for non-polar schemes
end

% figure('Name','Sampling Scheme'); scatter(sp(:,1),sp(:,2));

%% mask for the Fourier coefficients; default: full sampling 

M = length(sp);
mask = ones(M,1); 

%% initialize NUFFT 

disp('Initializing NUFFT...')
stx = nufft_init(2*pi*1/2^R*sp(:,1),2^R,7,2*2^R,0,'minmax:tuned');
sty = nufft_init(2*pi*1/2^R*sp(:,2),2^R,7,2*2^R,0,'minmax:tuned');
st = nufft_init(2*pi*1/2^R*sp,[2^R 2^R],[7 7],[2*2^R 2*2^R],[0 0],'minmax:tuned');

%% Precompute the Fourier transforms of the scaling function and the boundary (in both directions)
%  This step requires time

disp('Precomputing Fourier transforms of the scaling functions...')
if a==1  % the case of Haar wavelets
    ft_sca_x = haar_phi_ft(sp(:,1)/2^R); ft_sca_y = haar_phi_ft(sp(:,2)/2^R);
else     % the case of boundary corrected wavelets
    [ft_sca_L_x, ft_sca_x, ft_sca_R_x, ft_sca_L_y, ft_sca_y, ft_sca_R_y] = CDJV_Setup_NU2D(a, R, sp);
end

%% Calculate Fourier samples

disp('Calculating the Fourier samples...')

if K<=64 && strcmp(ss,'polar')
    omega_plain = Samples_NU(sp,'((1-x)*y+sin(y))^2', 0, 1, 0, 1);
    % %or:
    % omega_plain = Samples_NU(sp, '1', 0.25, 0.75, 0.25, 0.75 );
else
    [Xs,Ys] = meshgrid(0+0.5/6000:1/6000:1, 0+0.5/6000:1/6000:1);
    ffs = f(Xs,Ys).';
    S=length(ffs);
    sts = nufft_init(2*pi*1/S*sp,[S S],[10 10],[2*S 2*S],[0 0]);
    omega_plain = nufft((1/S)^2*ffs,sts);
end
omega = mask.*sqrt(mu).*omega_plain;

%% Compute the function GS_handle and solve for wavelet coefficient z

disp('Computing GS reconstruction...')
tGS_start = tic;  
if a==1
    GS_handle = @(x,mode) Haar_Op2_VecHandle_NU(x, mode, R, ft_sca_x, ft_sca_y, mu, st, mask(:));
else
    GS_handle = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty, mask(:));
end
z = lsqr(GS_handle,omega(:),10^(-12),50);
wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z))).';

%% reconstruction from wavelets

p=12;
if a==1
    reconim = get2DHaarReconstruction(wcoeff, p, R); 
else
    reconim = get2DReconstruction(wcoeff, a, p, R);
end

[x_supp,y_supp]=meshgrid(0:2^-p:1-2^-p);
orig=f(x_supp,y_supp);
figure('Name', 'Original function'); imshow(real(orig))
figure('Name', 'GS reconstruction'); imshow(real(reconim))

tGS = toc(tGS_start); 
tTotal = toc(tTotal_start); 
fprintf('Total computing time is %d minutes and %f seconds, out of which time for GS is %d minutes and %f seconds \n',floor(tTotal/60),rem(tTotal,60),floor(tGS/60),rem(tGS,60))

%% Calulate error of GS reconstruction and distance from the projection

error_gs=sqrt(sum(sum(abs(orig - reconim).^2))/2^(2*p));
fprintf('GS reconstruction error is %d \n',error_gs); 

title(['GS error: ', num2str(error_gs)]);

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
    w1 = wcoeff(2^j+1:2^(j+1), 2^j+1:2^(j+1));
    w2 = wcoeff(2^j+1:2^(j+1), 1:2^j);
    w3 = wcoeff(1:2^j, 2^j+1:2^(j+1)); 
    w_vec(ind(j-j0+1)+1:ind(j-j0+1)+ind(j-j0+2)) = [w1(:) ;w2(:) ; w3(:)];
end
figure('Name', 'Reconstructed wavelet coefficients'); stem(abs(w_vec)); set(gca,'xscal','log')

%% Compute the gridding reconstruction (truncated weighted Fourier series)

disp('Computing Gridding reconstruction...')
stg=nufft_init(2*pi*(1/2^p)*sp,[2^p 2^p],[7 7],[2*2^p 2*2^p],[0 0]);
f_gr=nufft_adj(1/4*mu.*omega_plain,stg).';

error_gr=sqrt(sum(sum(abs(orig - f_gr).^2))/2^(2*p));
fprintf('Error of Gridding reconstruction is %d \n', error_gr)

figure('Name', 'Gridding reconstruction'); imshow(real(f_gr))
title(['Gridding error: ', num2str(error_gr)]);

%%
clear('j','ind','w1','w2','w3')

% Copyright (c) 2014. Clarice Poon and Milana Gataric
