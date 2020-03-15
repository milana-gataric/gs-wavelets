%% 1d wavelet reconstruction when Fourier samples are taken on a jittered sampling scheme; 
%  default: full sampling & lsqr

clear all;
close all;

%% Parameters

jitter = 0.1;
epsilon = 0.75; 
delta = epsilon+2*jitter; % density of the sampling points; delta<1

a = 4; % number of vanishing moments \in {1,2,3,4,5,6,7,8}
R = 7; % maximum scale of wavelet coefficients

N = 2^R; % number of reconstruction vectors
K = N; % band-width of the samples, i.e. samples are taken from the interval [-K,K]

f =  @(x) cos(3*pi*x)+sin(6*pi*x).*(x+1/2).^2;   % function to be reconstructed
% % or:
% f = @(x) -exp(cos(6*pi*x)+sin(4*pi*x)).*cos(10*pi*x)+cos(4*pi*x);
% % or...

%% Sampling scheme and density compensation factors (Voronoi regions)

[sp,M] = jittered_sampling_1d(jitter,epsilon,K);
figure('Position', [200, 300, 1200, 670],'Name','1D reconstruction from jittered Fourier samples');
subplot(2,2,1); plot(sp,0); title('Sampling scheme');

fprintf('Using %d DB%d wavelets and %d jittered Fourier samples for GS reconstruction...\n',2^R,a,M);

mu=0.5*[2*(sp(2)-sp(1)); sp(3:M)-sp(1:M-2); 2*(sp(M)-sp(M-1))]; 

%% initialize NUFFT 

st = nufft_init(2*pi*1/N*sp,N,10,2*N,0,'minmax:tuned');

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

if a==1  % case of Haar wavelets
    ft_sca= haar_phi_ft(sp./2^R);
else  % case of boundary corrected wavelets
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup_NU1D(a, R, sp);
end
%% mask for the Fourier coefficients; default: full sampling 

mask = ones(M,1); 

%% Calculate Fourier samples

omega_plain = integral(@(x)f(x).*exp(-1i*2*pi*sp*x),0,1,'ArrayValued',true,'AbsTol',1e-12);
omega = mask.*sqrt(mu).*omega_plain;

%% Compute the Generalized function GS_handle and solve for wavelet coefficient z
if a==1
    GS_handle = @(x,mode) Haar_Op_Handle_NU(mode, x, R, ft_sca, mu, st, mask);
else
    GS_handle = @(x,mode) CDJV_Op_Handle_NU(mode, x, R, a, ft_sca, ft_sca_L, ft_sca_R, sp, mu, st, mask);
end
z = lsqr(GS_handle,omega(:),10^(-12));

%% reconstruction from wavelets

p=16;
if a==1
    recons = get1DHaarReconstruction(z, p, R);
else
    recons = get1DReconstruction(z, a, p, R);
end

orig=f((2^-p:2^-p:1)');
error_gs = sqrt(sum(abs(orig-recons).^2)/2^(p));
fprintf('GS error is %d \n',error_gs)

subplot(2,2,2); plot(2^-p:2^-p:1, orig); title('Original function'); axis([-0.02, 1.02 min(orig)-0.2 max(orig)+0.2])
subplot(2,2,3); plot(2^-p:2^-p:1, real(recons), 'm'); title({'GS reconstruction'; ['Error: ',num2str(error_gs)]}); axis([-0.02, 1.02 min(orig)-0.2 max(orig)+0.2])

%% Gridding reconstruction

stg=nufft_init(2*pi*(1/2^p)*sp,2^p,10,2*2^p,0); 
y=nufft_adj(mu.*omega_plain,stg);
 
error_gr = sqrt(sum(abs(f((2^-p:2^-p:1))-y.').^2)/2^(p));
fprintf('Gridding error is %d \n',error_gr)

subplot(2,2,4); plot(2^-p:2^-p:1, real(y)); title({'Gridding reconstruction'; ['Error: ',num2str(error_gr)]}); axis([-0.02, 1.02 min(orig)-0.2 max(orig)+0.2])


% Copyright (c) 2014. Clarice Poon and Milana Gataric
