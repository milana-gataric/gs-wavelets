%% 1d reconstruction matrix for Haar scaling functions when Fourier samples are 
%  taken on a nonuniform sampling scheme: radial (polar) sampling scheme or spiral;

clear all;

%% Parameters

D = sqrt(2)/4; % the bound for the density in the Euclidean norm 
               % i.e. delta < D; max D is sqrt(2)/4
ss = 'polar'; % type of the sampling scheme, it can be a polar or an (interleaving) 
              % spiral: 'polar', 'spiral 1', 'spiral 2', 'spiral 3', ...
               
R = 5; % maximum scale of wavelet coefficients

N = 2^R; % number of reconstruction vectors
K = N; % band-width of the samples

%% Sampling scheme and density compensation factors 

disp('Calculating the sampling scheme and density compensation factors...')

if exist('mri_density_comp', 'file') % the same function has different names across the different versions of the package IR
    density_fn = @(x) mri_density_comp(x, 'voronoi');
else
    density_fn = @(x) ir_mri_density_comp(x, 'voronoi');
end

if strcmp(ss,'polar')
    sp = polar_sampling(K, D);
    mu = 4*density_fn(sp); % takes time for large reconstructions
elseif strfind(ss, 'spiral')
    [~, leaves] = strtok(ss);
    sp = spiral_sampling(K, D, str2double(leaves));
    mu = 4*density_fn(sp); % takes time for large reconstructions
    mu(mu>5)=5; % a correction is needed since the function 'mri_density_comp' doesn't work perfectly for non-polar schemes
end

%% initialize NUFFT 
% 
% disp('Initializing NUFFT...')
% stx = nufft_init(2*pi*1/2^R*sp(:,1),2^R,7,2*2^R,0,'minmax:tuned');
% sty = nufft_init(2*pi*1/2^R*sp(:,2),2^R,7,2*2^R,0,'minmax:tuned');
% st = nufft_init(2*pi*1/2^R*sp,[2^R 2^R],[7 7],[2*2^R 2*2^R],[0 0],'minmax:tuned');
% 
%% Precompute the fourier transforms of the scaling functions
% 
% disp('Precomputing Fourier transforms of the scaling functions...')
% ft_sca_x = haar_phi_ft(sp(:,1)/2^R); ft_sca_y = haar_phi_ft(sp(:,2)/2^R);

%% Computation of the GS matrix

disp('Cmputing the GS matrix...');
diag_haar=sqrt(mu).*1/N.*exp(-1i*pi*1/N.*(sp(:,1)+sp(:,2))).*sinc(1/N.*sp(:,1)).*sinc(1/N.*sp(:,2));
m=0; U=zeros(length(sp),N^2);
for m1=0:N-1
    for m2=0:N-1
        m=m+1;
        U(:,m)=diag_haar(:).*exp(-1i*2*pi*1/N*(m1*sp(:,1)+m2*sp(:,2)));
    end
end

min_sing = sqrt(min(eig(U'*U)));
disp('Minimal singular eigenvalue is')
disp(min_sing)

rec_const = (1+1)/min_sing;
disp('Reconstruction constant estimate is')
disp(rec_const)
% Copyright (c) 2014. Clarice Poon and Milana Gataric
