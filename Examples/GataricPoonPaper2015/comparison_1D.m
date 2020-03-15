% 1d wavelet reconstruction and truncated Fourier of a continuous function 
% for different sampling schemes: uniform, jittered and log

% comment: shows superior reconstruction by GS for different sampling shemes

clear all;
close all;

%% Parameters

a = 4; % number of vanishing moments
R = 6; % maximum scale of wavelet coefficients

N = 2^R; % number of reconstruction vectors
K = N; % band-width of the samples, i.e. samples are taken from the interval [-K,K]

f = @(x) -exp(x.*cos(4*pi*x)).*cos(7*pi*x)+sin(3*pi*x); % function to be reconstructed
xylims=[-0.02 1.02 -2.6 3.4];

%f = @(x) x.^2.*cos(7*pi*x); % function to be reconstructed
%xylims=[-0.02 1.02 -1.2 1];

p = 18; % this resolution takes more time

figure('Name','Original function'); plot(2^-p:2^-p:1,f(2^-p:2^-p:1)); axis(xylims)

figure('Position', [50, 100, 1600, 670],'Name','1D smooth function recovery from different sampling schemes');

%% Uniform sampling

disp('Uniform sampling...')

eps = 1;
S=1;
M = S*2^R/eps;
sp=((M:-1:-M+1)*eps)';

fprintf('Number of samples used is %d \n',length(sp))

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

W = 1/(2^R)*sp;
[ft_sca_L, ft_sca, ft_sca_R] = CDJV_FT_Vec(a, W);

%% Fourier samples

omega_plain = integral(@(x)f(x).*exp(-1i*2*pi*sp*x),0,1,'ArrayValued',true,'AbsTol',1e-14);
omega = omega_plain;

%% Compute the Generalized function GS_handle and solve for wavelet coefficient z

GS_handle = @(x,mode) CDJV_Op_Handle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
z = lsqr(GS_handle,omega(1:2*M),10^(-12));

%% reconstruction from wavelets

recons = get1DReconstruction(z, a, p, R);

error_gs_un = sqrt(sum(abs(f((0:2^-p:1-2^-p)')-recons).^2)/2^(p));

fprintf('GS error from uniform samples is %d \n',error_gs_un)

subplot(2,3,1); plot(2^-p:2^-p:1, real(recons),'m'); title({'GS reconstruction from uniform samples'; ['Error: ',num2str(error_gs_un)]}); axis(xylims)

%% Truncated Fourier representation

Supp=2^-p:2^-p:1;
y=zeros(size(Supp));
for j=M:-1:-M+1
    y=y+sqrt(eps)*exp(2*eps*1i*j*Supp*pi)*omega(M-j+1);
end

error_ft_un = sqrt(sum(abs(f((2^-p:2^-p:1))-y).^2)/2^(p));

fprintf('Truncated Fourier error from uniform samples is %d \n',error_ft_un)

subplot(2,3,4); plot(Supp, real(y)); title({'Truncated Fourier reconstruction from uniform samples'; ['Error: ',num2str(error_ft_un)]}); axis(xylims)   

%% Jittered samples

disp('Jittered sampling...')

jitter = 0.1;
epsilon = 0.77; 
delta = epsilon+2*jitter; % density of the sampling points; delta<1

[sp,M] = jittered_sampling_1d(jitter,epsilon,K);

mu=0.5*[2*(sp(2)-sp(1)); sp(3:M)-sp(1:M-2); 2*(sp(M)-sp(M-1))]; 

fprintf('Number of samples used is %d \n',length(sp))

%% initialize NUFFT 

st = nufft_init(2*pi*1/N*sp,N,10,2*N,0,'minmax:tuned');

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

W = 1/(2^R)*sp;
[ft_sca_L, ft_sca, ft_sca_R] = CDJV_FT_Vec(a, W);

%% Define a function to be reconstructed and calculate its Fourier samples

omega_plain = integral(@(x)f(x).*exp(-1i*2*pi*sp*x),0,1,'ArrayValued',true,'AbsTol',1e-14);
omega = sqrt(mu).*omega_plain;

%% Compute the Generalized function GS_handle and solve for wavelet coefficient z

GS_handle = @(x,mode) CDJV_Op_Handle_NU(mode, x, R, a, ft_sca, ft_sca_L, ft_sca_R, sp, mu, st);
z = lsqr(GS_handle,omega(:),10^(-12));

%% reconstruction from wavelets

recons = get1DReconstruction(z, a, p, R);

error_gs_jit = sqrt(sum(abs(f((2^-p:2^-p:1)')-recons).^2)/2^(p));
fprintf('GS error from jittered samples is %d \n',error_gs_jit)

subplot(2,3,2); plot(2^-p:2^-p:1, real(recons), 'm'); title({'GS reconstruction from jittered samples'; ['Error: ',num2str(error_gs_jit)]}); axis(xylims)


%% Gridding reconstruction

stg=nufft_init(2*pi*(1/2^p)*sp,2^p,10,2*2^p,0); 
y=nufft_adj(mu.*omega_plain,stg);

error_gr_jit= sqrt(sum(abs(f((2^-p:2^-p:1))-y.').^2)/2^(p));
fprintf('Gridding error from jittered samples is %d \n',error_gr_jit)

subplot(2,3,5); plot(2^-p:2^-p:1, real(y)); title({'Gridding reconstruction from jittered samples'; ['Error: ',num2str(error_gr_jit)]}); axis(xylims)

%% Log sampling

disp('Log sampling...')

delta = 0.97; % density of the sampling points; delta<1
[sp,M]=log_sampling(delta,K);

mu=0.5*[2*(sp(2)-sp(1)); sp(3:M)-sp(1:M-2); 2*(sp(M)-sp(M-1))]; 

fprintf('Number of samples used is %d \n',length(sp))

%% initialize NUFFT 

st = nufft_init(2*pi*1/N*sp,N,10,2*N,0,'minmax:tuned');

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

W = 1/(2^R)*sp;
[ft_sca_L, ft_sca, ft_sca_R] = CDJV_FT_Vec(a, W);


%% calculate its Fourier samples

omega_plain = integral(@(x)f(x).*exp(-1i*2*pi*sp*x),0,1,'ArrayValued',true,'AbsTol',1e-14);
omega = sqrt(mu).*omega_plain;

%% Compute the Generalized function GS_handle and solve for wavelet coefficient z

GS_handle = @(x,mode) CDJV_Op_Handle_NU(mode, x, R, a, ft_sca, ft_sca_L, ft_sca_R, sp, mu, st);
z = lsqr(GS_handle,omega(:),10^(-12));

%% reconstruction from wavelets

recons = get1DReconstruction(z, a, p, R);

error_gs_log = sqrt(sum(abs(f((2^-p:2^-p:1)')-recons).^2)/2^(p));
fprintf('GS error from log samples is %d \n',error_gs_log)

subplot(2,3,3); plot(2^-p:2^-p:1, real(recons), 'm'); title({'GS reconstruction from log samples'; ['Error: ',num2str(error_gs_log)]}); axis(xylims)


%% Gridding reconstruction

stg=nufft_init(2*pi*(1/2^p)*sp,2^p,10,2*2^p,0); 
y=nufft_adj(mu.*omega_plain,stg);

error_gr_log = sqrt(sum(abs(f((2^-p:2^-p:1))-y.').^2)/2^(p));
fprintf('Gridding error from log samples is %d \n',error_gr_log)

subplot(2,3,6); plot(2^-p:2^-p:1, real(y)); title({'Gridding reconstruction from log samples'; ['Error: ',num2str(error_gr_log)]}); axis(xylims)
 





% Copyright (c) 2014. Clarice Poon and Milana Gataric
