%% Generates 2d boundary corrected wavelet matrix when Fourier samples are 
%  taken on a nonuniform sampling scheme: radial (polar) sampling scheme or spiral;
%  Also, an estamation of the reconstruction constant is calculated

% Warning: Precomputation steps take longer time for K>=64 (for K=64, 
% about 10 minutes); use precomputed data for bigger examples

clear all;

%% Parameters

D = sqrt(2)/4; % the bound for the density in the Euclidean norm 
               % i.e. delta < D; max D is sqrt(2)/4
ss = 'polar'; % type of the sampling scheme, it can be a polar or an (interleaving) 
              % spiral: 'polar', 'spiral 1', 'spiral 2', 'spiral 3', ...
               
a = 2; % number of vanishing moments
R = 4; % maximum scale of wavelet coefficients

N = 2^R; % number of reconstruction vectors
K = N; % band-width of the samples

%% Sampling scheme (sampling points, number of sampling points) and 
%  density compensation factors (Voronoi regions)

disp('Calculating the sampling scheme and density compensation factors...')

if exist('mri_density_comp', 'file') % the same function has different names across the different versions of the package IR
    density_fn = @(x) mri_density_comp(x, 'voronoi');
else
    density_fn = @(x) ir_mri_density_comp(x, 'voronoi');
end

if strcmp(ss,'polar')
    sp = polar_sampling(K, D);
    mu = 4*density_fn(sp);
elseif strfind(ss, 'spiral')
    [~, leaves] = strtok(ss);
    sp = spiral_sampling(K, D, str2double(leaves));
    mu = 4*density_fn(sp); 
    mu(mu>5)=5; % a correction is needed since the function 'mri_density_comp' doesn't work perfectly for non-polar schemes
end

%% initialize NUFFT 

disp('Initializing NUFFT...')
stx = nufft_init(2*pi*1/2^R*sp(:,1),2^R,7,2*2^R,0,'minmax:tuned');
sty = nufft_init(2*pi*1/2^R*sp(:,2),2^R,7,2*2^R,0,'minmax:tuned');
st = nufft_init(2*pi*1/2^R*sp,[2^R 2^R],[7 7],[2*2^R 2*2^R],[0 0],'minmax:tuned');

%% The fourier transforms of the scaling functions

disp('Precomputing Fourier transforms of the scaling functions...')
[ft_sca_L_x, ft_sca_x, ft_sca_R_x, ft_sca_L_y, ft_sca_y, ft_sca_R_y] = CDJV_Setup_NU2D(a, R, sp);

%% Computation of the GS matrix

disp('Computing the GS matrix...')

U=zeros(length(sp),2^(2*R));
wc= zeros(2^R,2^R); 

index = 1;

for u=1:2^R
    for v=1:u

        wc(u,v)=1;
        y  = CDJV_Op2_Handle_NU(1, wc, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);
        U(:, index) = y(:);
        index=index+1;
        wc(u,v)=0;
        
        if u~= v
            wc(v,u)=1;
            y  = CDJV_Op2_Handle_NU(1, wc, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);

    
            U(:, index) = y(:);
            index=index+1;
            wc(v,u)=0;
        end

       
    end
end

min_sing = sqrt(min(eig(U'*U)));
disp('Minimal singular eigenvalue is')
disp(min_sing)

rec_const = (1+1)/min_sing;
disp('Reconstruction constant estimate is')
disp(rec_const)
% Copyright (c) 2014. Clarice Poon and Milana Gataric
