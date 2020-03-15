% illustration of error decrease for 1d reconstruction with different wavelets 
% when Fourier samples are taken uniformly 

clear all;
close all;

eps = 1; % distance between consecutive samples; require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; % Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1

f = @(x) x.*cos(3*pi*x); % function to be reconstructed
xylims=[-0.02 1.02 -1.2 1];

p = 18;

figure('Name','Original function'); plot(2^-p:2^-p:1,f(2^-p:2^-p:1)); axis(xylims)

Rrange=5:8;
R0=Rrange(1);

%% Haar

a = 1; % number of vanishing moments


for R=Rrange % maximum scale of wavelet coefficients
    
    fprintf('%d Haar: \n',2^R);
    M = S*2^R/eps;
    
    omega = integral(@(x)f(x).*exp(-1i*2*pi*((M:-1:-M+1)*eps)'*x),0,1,'ArrayValued',true,'AbsTol',1e-14);
    
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');

    GS_handle = @(x,mode) Haar_Op_Handle(mode, x, S, eps, R, ft_sca);
    z = lsqr(GS_handle,omega(1:2*M),10^(-12));
    
    recons = get1DHaarReconstruction(z, p, R);

    error_gs_Haar(R-R0+1) = sqrt(sum(abs(f((0:2^-p:1-2^-p)')-recons).^2)/2^(p)); %#ok<*SAGROW>
    log_Haar_err(R-R0+1) = -log2(error_gs_Haar(R-R0+1))/log2(2*M);
    
    fprintf('GS error is %d \n',error_gs_Haar(R-R0+1));
    fprintf('error/#samples^(-a) is %d \n',error_gs_Haar(R-R0+1)/(2*M)^(-a));
    fprintf('-log2(error)/log2(#samples) is %d \n',log_Haar_err(R-R0+1));
    disp(' ')
    
       
end

%% DB2, DB3

for a=2:3 % number of vanishing moments
    for R=Rrange % maximum scale of wavelet coefficients
        fprintf('%d DB%d: \n',2^R, a)
        
        M = S*2^R/eps;

        omega = integral(@(x)f(x).*exp(-1i*2*pi*((M:-1:-M+1)*eps)'*x),0,1,'ArrayValued',true,'AbsTol',1e-14);

        [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);

        GS_handle = @(x,mode) CDJV_Op_Handle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
        z = lsqr(GS_handle,omega(1:2*M),10^(-12));
        
        recons = get1DReconstruction(z, a, p, R, 24); 
        
        error_gs_DB(R-R0+1,a-1) = sqrt(sum(abs(f((0:2^-p:1-2^-p)')-recons).^2)/2^(p)); %#ok<*SAGROW>
        log_DB_err(R-R0+1,a-1)  = -log2(error_gs_DB(R-R0+1,a-1))/log2(2*M);

        fprintf('GS error is %d \n', error_gs_DB(R-R0+1,a-1));
        fprintf('error/#samples^(-a) is %d \n',error_gs_DB(R-R0+1,a-1)/(2*M)^(-a));
        fprintf('-log2(error)/log2(#samples) is %d \n',log_DB_err(R-R0+1,a-1));
        disp(' ')
    end
end

for R=Rrange 
    
    M = 2^R;
    fprintf('Truncated Fourier from %d points: \n',2*M);
        
    omega = integral(@(x)f(x).*exp(-1i*2*pi*((M:-1:-M+1)*eps)'*x),0,1,'ArrayValued',true,'AbsTol',1e-14);
    
    Supp=2^-p:2^-p:1;
    y=zeros(size(Supp));
    for j=M:-1:-M+1
        y=y+sqrt(eps)*exp(2*eps*1i*j*Supp*pi)*omega(M-j+1);
    end
    
    error_ft(R-R0+1) = sqrt(sum(abs(f((2^-p:2^-p:1))-y).^2)/2^(p));
    log_ft(R-R0+1) = -log2(error_ft(R-R0+1))/log2(2*M);
    
    fprintf('TF error is %d \n', error_ft(R-R0+1));
    fprintf('error/#samples is %d \n',error_ft(R-R0+1)/(2*M));
    fprintf('-log2(error)/log2(#samples) is %d \n',log_ft(R-R0+1));
    disp(' ')
end



figure('Name','-log2(error)');
plot(Rrange+1,-log2(error_ft),'ok-','LineWidth',1,'MarkerSize',5);
hold on;
plot(Rrange+1,-log2(error_gs_Haar),'sb-','LineWidth',1,'MarkerSize',5);
plot(Rrange+1,-log2(error_gs_DB(:,1)),'dm-','LineWidth',1,'MarkerSize',5);
plot(Rrange+1,-log2(error_gs_DB(:,2)),'^g-','LineWidth',1,'MarkerSize',5);
set(gca,'XTick',Rrange+1)
xlabel('log2(#samples)');
ylabel('-log2(error)');
l= legend('TFS','Haar','DB2','DB3');
set(l,'Interpreter','Latex','FontSize',13, 'Location','NorthWest');





% Copyright (c) 2014. Clarice Poon and Milana Gataric
