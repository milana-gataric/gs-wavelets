%% 1d reconstruction matrix for Haar scaling functions and a jittered
%  sampling scheme

clear all;

%% Parameters

jitter = 0.1;
epsilon = 0.75; 
delta = epsilon+2*jitter; % density of the sampling points; delta<1

R = 8; % maximum scale of wavelet coefficients

N = 2^R; % number of reconstruction vectors
K = N; % band-width of the samples, i.e. samples are taken from the interval [-K,K]

%% Sampling scheme (sampling points, number of sampling points) and 
%  density compensation factors (Voronoi regions)

[sp,M] = jittered_sampling_1d(jitter,epsilon,K);

mu=0.5*[2*(sp(2)-sp(1)); sp(3:M)-sp(1:M-2); 2*(sp(M)-sp(M-1))]; 

%% initialize NUFFT 

% st = nufft_init(2*pi*1/N*sp,N,10,2*N,0,'minmax:tuned');
 
%% The fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

% ft_sca= haar_phi_ft(sp./2^R);

%% Computation of the matrix

U=zeros(M,N);

n=0:N-1;

for k=1:M
    U(k,:)=sqrt(mu(k))*(1/sqrt(2)).*(sqrt(2/N)).*exp(-1i*pi*sp(k)*1/N).*sinc(sp(k)*1/N).*exp(-1i*pi*sp(k).*(2.*n/N-1)); 
end

min_sing = sqrt(min(eig(U'*U)));
disp('Minimal singular eigenvalue is')
disp(min_sing)

rec_const = (1+delta)/min_sing;
disp('Reconstruction constant estimate is')
disp(rec_const)

% Copyright (c) 2014. Clarice Poon and Milana Gataric
