close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%
%%%Parameters
%%%%%%%%%%%%%%%%%%%%%%
a=8; %number of vanishing moments
R = 5;
eps = 1/2; %require 1/eps to be a natural number
M=2^R/eps;

%%Define a function on [0,1] to be reconstructed, and the Fourier samples
%%up to frequency M
f = @(x) x.^3+1;
fn_string = 'x^3+1';
omega = Samples( eps, M, fn_string, 0, 1);




%%GS l1 reconstruction
S=1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R_ext= max(R+2,10); %maximum scale of wavelet coefficients
M_ext = 2^R_ext/eps;
mask = zeros(2*M_ext,1);
mask(M_ext+1-M:M_ext+M)=1;
omega_ext = zeros(2*M_ext,1);
omega_ext(M_ext+1-M:M_ext+M)=omega;

if a==1 %%%For HAAR
    ft_sca_ext = haar_phi_ft(((M_ext:-1:-M_ext+1)./(2^R_ext/eps))');
    GS_l1Op = @(x,mode) Haar_Op_Handle(mode, x, S, eps, R_ext, ft_sca_ext, mask);
else %%%For CDJV
    [ft_sca_L_ext, ft_sca_ext, ft_sca_R_ext] = CDJV_Setup(a, eps, S, R_ext);
    GS_l1Op = @(x,mode) CDJV_Op_Handle(mode, x, S, eps, R_ext, a, ft_sca_ext, ft_sca_L_ext, ft_sca_R_ext, mask);
end
disp('Computing wavelet reconstruction using GS l1...');
opts = spgSetParms('verbosity',1, 'bpTol', 1e-10,'optTol', 1e-10, 'iterations', 1000);
wcoeff_consistent = spg_bp(GS_l1Op,omega_ext,opts);



%% GS LSQR reconstruction

if a==1
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');
    GS_lsqrOp = @(x,mode) Haar_Op_Handle(mode, x, S, eps, R, ft_sca);
else
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);
    GS_lsqrOp = @(x,mode) CDJV_Op_Handle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
end
wcoeff_lsqr = lsqr(GS_lsqrOp,omega, 1e-10);




%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
p=max(14, R+1);
if a==1
    
    recons_gsl1 = get1DHaarReconstruction(wcoeff_consistent, p, R_ext);
    recons_lsqr = get1DHaarReconstruction(wcoeff_lsqr, p, R);
    
    
else
    recons_gsl1 = get1DReconstruction(wcoeff_consistent, a, p, R_ext);
    recons_lsqr = get1DReconstruction(wcoeff_lsqr, a, p, R);
end

x_supp=(0:2^-p:1-2^-p)';

orig=f(x_supp);
figure('Position', [150, 150, 1250, 500],'Name','1D reconstruction of a discontinuous funcion from sparse uniform samples');
subplot(1,4,1); plot(x_supp,real(orig)); title('Original function');
subplot(1,4,2); plot(x_supp,real(recons_gsl1)); title('GS l1 reconstruction');
subplot(1,4,3); plot(x_supp,real(recons_lsqr)); title('GS lsqr reconstruction');
gsl1_error = norm(f(x_supp)-recons_gsl1)/2^p;
lsqr_error = norm(f(x_supp)-recons_lsqr)/2^p;





%% Truncated Fourier representation
fourier_rep=zeros(size(x_supp));
for j=M:-1:-M+1
    fourier_rep=fourier_rep+sqrt(eps)*exp(2*eps*1i*j*x_supp*pi)*omega(M-j+1);
end

subplot(1,4,4); plot(x_supp,real(fourier_rep)); title('Truncated Fourier reconstruction');

%
fourier_error = sqrt(sum(abs(f(x_supp)-fourier_rep).^2)/2^(p));

fprintf('Consistency of GS l1 with given samples is %d\n', norm(omega_ext - GS_l1Op(wcoeff_consistent,1)));

fprintf('Consistency of GS lsqr with given samples is %d\n', norm(omega - GS_lsqrOp(wcoeff_lsqr,1)));
fprintf('Error from GS l1 reconstruction is %d\n', gsl1_error);
fprintf('Error from GS lsqr reconstruction is %d\n', lsqr_error);
fprintf('Error from truncated Fourier reconstruction is %d\n', fourier_error);


% Copyright (c) 2014. Clarice Poon and Milana Gataric
