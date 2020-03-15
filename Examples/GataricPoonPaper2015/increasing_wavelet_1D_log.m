% illustration of error decrease for 1d reconstruction with different wavelets 
% when Fourier samples are taken on log sampling scheme

clear all;
close all;

delta = 0.97; % density of the sampling points; delta<1

f = @(x) x.*cos(3*pi*x); % function to be reconstructed
xylims=[-0.02 1.02 -1.2 1];

p = 18;

figure('Name','Original function'); plot(2^-p:2^-p:1,f(2^-p:2^-p:1)); axis(xylims)

Rrange=[5 6 7 8];
R0=Rrange(1);

%% Haar

a = 1; % number of vanishing moments

for R=Rrange % maximum scale of wavelet coefficients
    fprintf('%d Haar: \n',2^R);
    
    N = 2^R;
    K = N;
    
    [sp,M]=log_sampling(delta,K);
    mu=0.5*[2*(sp(2)-sp(1)); sp(3:M)-sp(1:M-2); 2*(sp(M)-sp(M-1))]; 
    
    omega = sqrt(mu).*integral(@(x)f(x).*exp(-1i*2*pi*sp*x),0,1,'ArrayValued',true,'AbsTol',1e-14);
    
    ft_sca= haar_phi_ft(sp./2^R);
    
    st = nufft_init(2*pi*1/N*sp,N,10,2*N,0,'minmax:tuned');

    GS_handle = @(x,mode) Haar_Op_Handle_NU(mode, x, R, ft_sca, mu, st);
    z = lsqr(GS_handle,omega(:),10^(-12));
    
    recons = get1DHaarReconstruction(z, p, R);

    error_gs_Haar(R-R0+1) = sqrt(sum(abs(f((0:2^-p:1-2^-p)')-recons).^2)/2^(p)); %#ok<*SAGROW>

    fprintf('GS error is %d \n',error_gs_Haar(R-R0+1));
    fprintf('error/K^(-a) is %d \n',error_gs_Haar(R-R0+1)/(K)^(-a));
    fprintf('-log2(error)/log2(K) is %d \n',-log2(error_gs_Haar(R-R0+1))/log2(K));
    disp(' ')
end

%% DB2, DB3

for a=2:3 % number of vanishing moments
    for R=Rrange % maximum scale of wavelet coefficients
        fprintf('%d DB%d: \n',2^R, a)
        
        N = 2^R;
        K = N;
    
        [sp,M]=log_sampling(delta,K);
        mu=0.5*[2*(sp(2)-sp(1)); sp(3:M)-sp(1:M-2); 2*(sp(M)-sp(M-1))]; 
        
        omega = sqrt(mu).*integral(@(x)f(x).*exp(-1i*2*pi*sp*x),0,1,'ArrayValued',true,'AbsTol',1e-14);

        W = 1/(2^R)*sp;
        [ft_sca_L, ft_sca, ft_sca_R] = CDJV_FT_Vec(a, W);
        
        st = nufft_init(2*pi*1/N*sp,N,10,2*N,0,'minmax:tuned');

        GS_handle = @(x,mode) CDJV_Op_Handle_NU(mode, x, R, a, ft_sca, ft_sca_L, ft_sca_R, sp, mu, st);
        z = lsqr(GS_handle,omega(:),10^(-12),30);
        
        recons = get1DReconstruction(z, a, p, R,24);

        error_gs_DB(R-R0+1,a-1) = sqrt(sum(abs(f((0:2^-p:1-2^-p)')-recons).^2)/2^(p)); %#ok<*SAGROW>
        
        fprintf('GS error is %d \n', error_gs_DB(R-R0+1,a-1));
        fprintf('error/K^(-a) is %d \n',error_gs_DB(R-R0+1,a-1)/(K)^(-a));
        fprintf('-log2(error)/log2(K) is %d \n',-log2(error_gs_DB(R-R0+1,a-1))/log2(K));
        disp(' ')
    end
end

%% gridding

for R=Rrange 
    
    K = 2^R;
    fprintf('Gridding from bandwidth %d: \n',K);
    
    [sp,M]=log_sampling(delta,K);
    mu=0.5*[2*(sp(2)-sp(1)); sp(3:M)-sp(1:M-2); 2*(sp(M)-sp(M-1))]; 
    
    omega = integral(@(x)f(x).*exp(-1i*2*pi*sp*x),0,1,'ArrayValued',true,'AbsTol',1e-14);
    
    stg=nufft_init(2*pi*(1/2^p)*sp,2^p,10,2*2^p,0); 
    y=nufft_adj(mu.*omega,stg);

    error_gr(R-R0+1) =sqrt(sum(abs(f((2^-p:2^-p:1))-y.').^2)/2^(p));
    log_gr(R-R0+1) = -log2(error_gr(R-R0+1))/log2(K);
    
    fprintf('Gridding error is %d \n', error_gr(R-R0+1));
    fprintf('error/K is %d \n',error_gr(R-R0+1)/(K));
    fprintf('-log2(error)/log2(K) is %d \n',log_gr(R-R0+1));
    disp(' ')
end

%% plots 

figure('Name','-log2(error) : log sampling');
plot(Rrange,-log2(error_gr),'ok-','LineWidth',1,'MarkerSize',5);
hold on;
plot(Rrange,-log2(error_gs_Haar),'sb-','LineWidth',1,'MarkerSize',5);
plot(Rrange,-log2(error_gs_DB(:,1)),'dm-','LineWidth',1,'MarkerSize',5);
plot(Rrange,-log2(error_gs_DB(:,2)),'^g-','LineWidth',1,'MarkerSize',5);
set(gca,'XTick',Rrange)
xlabel('log2(K)');
ylabel('-log2(error)');
l= legend('Gr','Haar','DB2','DB3');
set(l,'Interpreter','Latex','FontSize',13, 'Location','NorthWest');





% Copyright (c) 2014. Clarice Poon and Milana Gataric
