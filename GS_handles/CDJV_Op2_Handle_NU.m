% Applies 2D GS operator with boundary wavelets to a matrix when samples 
% are taken nonuniformly a form of a handle function y = G(mode,wc) 
%
% Output:  
%       y = G(mode,x) is  y = G*wc           if mode == 1 (or 'notransp');
%                         y = adjoint(G)*wc  if mode == 2 (or 'transp')
%       where G is the generalized sampling operator 
%       G : wc [a matrix of 2^R times 2^R of 2D wavelet coeff] -> y [a vector of M nonuniform Fourier samples] 
%
% Input:
%       mode: argument of the function handle; 
%             mode == 1 (or 'notransp') computes the forward operation and 
%             mode == 2 (or 'transp') computes the adjoint operation
%       x:    argument of the handle function;
%             a vector of (2^R)^2 2D wavelet coefficients for mode == 1 and
%             a vector of length(sp) NONUNIFORM Fourier coefficients for mode == 2
%       R: maximum scale of wavelet coefficients
%       a: number of vanishing moments
%       ft_sca_x, ft_sca_L_x, ft_sca_R_x: Fourier transform of the main, left, 
%                                         and right scaling function at sp(:,1)
%       ft_sca_y, ft_sca_L_y, ft_sca_R_y: Fourier transform of the main, left, 
%                                         and right scaling function at sp(:,2)
%       sp: a M times 2 matrix of sampling points
%       mu: a vector of M density compensetaion factors
%       st, stx, sty: initialize structure for NUFFT
%       mask: an optional argument which represents a mask for the Fourier coefficients

function y_concat = CDJV_Op2_Handle_NU(mode, wc, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, SamplePoints, mu, st, stx, sty, mask)

    L=2^R;
    
    if or(strcmp(mode,'notransp'),mode==1)
       
        x = IWT2_CDJV_noP(real(wc),ceil(log2(2*a))+1,a) +1i * IWT2_CDJV_noP(imag(wc),ceil(log2(2*a))+1,a);
        x_pad = zeros(L, L);
        x_pad(2:2^R-2*a+1, 2:2^R-2*a+1) = x(a+1:2^R-a, a+1:2^R-a);
        x_pad=x_pad.';
        
        y_concat = nufft(x_pad,st);
        
        y_concat = ft_sca_x.*(ft_sca_y.*y_concat);
        
        LxLy = zeros(length(SamplePoints),1); LxMy = zeros(length(SamplePoints),1); LxRy = zeros(length(SamplePoints),1);
        MxLy = zeros(length(SamplePoints),1); MxRy = zeros(length(SamplePoints),1);
        RxLy = zeros(length(SamplePoints),1); RxMy = zeros(length(SamplePoints),1); RxRy = zeros(length(SamplePoints),1);
        
        VL=zeros(2^R,a); VR=zeros(2^R,a); HU=zeros(a,2^R); HD=zeros(a,2^R); 
        VL(2:2^R-2*a+1, 1:a)=x(a+1:2^R-a, 1:a);
        VR(2:2^R-2*a+1, 1:a)=x(a+1:2^R-a, 2^R-a+1:2^R);
        HU(1:a, 2:2^R-2*a+1)=x(1:a,a+1:2^R-a);
        HD(1:a, 2:2^R-2*a+1)=x(2^R-a+1:2^R,a+1:2^R-a);
        for n=1:a
            for m=1:a
                LxLy = LxLy + ft_sca_L_y(:,m).*ft_sca_L_x(:,n)*x(m,n);
                LxRy = LxRy + ft_sca_R_y(:,m).*ft_sca_L_x(:,n).*exp(-1i*2*pi*SamplePoints(:,2))*x(2^R-a+m,n);
                RxLy = RxLy + ft_sca_L_y(:,m).*ft_sca_R_x(:,n).*exp(-1i*2*pi*SamplePoints(:,1))*x(m,2^R-a+n);
                RxRy = RxRy + exp(-1i*2*pi*(SamplePoints(:,1)+SamplePoints(:,2))).*ft_sca_R_y(:,m).*ft_sca_R_x(:,n)*x(2^R-a+m,2^R-a+n);
            end
            LxMy = LxMy + ft_sca_y.*ft_sca_L_x(:,n).*nufft(VL(:,n),sty);
            MxLy = MxLy + ft_sca_L_y(:,n).*ft_sca_x.*nufft(HU(n,:).',stx);
            MxRy = MxRy + exp(-1i*2*pi*SamplePoints(:,2)).*ft_sca_R_y(:,n).*ft_sca_x.*nufft(HD(n,:).',stx);
            RxMy = RxMy + ft_sca_y.*ft_sca_R_x(:,n).*exp(-1i*2*pi*SamplePoints(:,1)).*nufft(VR(:,n),sty);
        end
     
                         
        y_concat = y_concat + LxLy + LxMy + LxRy + MxLy + MxRy + RxLy + RxMy + RxRy;
        y_concat = y_concat.*sqrt(mu)*1/L;
        
        if nargin == 16
            y_concat = mask.*y_concat;
        end

    else
        x_trunc=zeros(2^R,2^R);
        wc=wc./L;
        wc = sqrt(mu).*wc;
        if nargin == 16
            wc = mask.*wc;
        end
        
        x = conj(ft_sca_y).*(conj(ft_sca_x).*wc);        
        xmm = nufft_adj(x,st);
        x_trunc(a+1:2^R-a,a+1:2^R-a) = xmm(2:2^R-2*a+1, 2:2^R-2*a+1).';
        
        for m=1:a
            for n=1:a
                x_trunc(m,n) = sum(conj(ft_sca_L_y(:,m)).*conj(ft_sca_L_x(:,n)).*wc);
                x_trunc(2^R-a+m, n) = sum(conj(ft_sca_R_y(:,m)).*conj(ft_sca_L_x(:,n)).*wc.*exp(1i*2*pi*SamplePoints(:,2)));
                x_trunc(m, 2^R-a+n) = sum(conj(ft_sca_L_y(:,m)).*conj(ft_sca_R_x(:,n)).*wc.*exp(1i*2*pi*SamplePoints(:,1)));
                x_trunc(2^R-a+m, 2^R-a+n) = sum(conj(ft_sca_R_y(:,m)).*conj(ft_sca_R_x(:,n)).*wc.*exp(1i*2*pi*(SamplePoints(:,1)+SamplePoints(:,2))));
            end
        end
        
        LxMy = zeros(2^R, a);  RxMy = zeros(2^R, a);
        MxLy = zeros(a,2^R); MxRy = zeros(a,2^R);
        for n=1:a 
            LxMy(:,n) = nufft_adj(conj(ft_sca_y).*conj(ft_sca_L_x(:,n)).*wc,sty);
            RxMy(:,n) = nufft_adj(conj(ft_sca_y).*conj(ft_sca_R_x(:,n)).*exp(1i*2*pi*(SamplePoints(:,1))).*wc,sty);
            MxLy(n,:) = nufft_adj(conj(ft_sca_L_y(:,n)).*conj(ft_sca_x).*wc,stx).';
            MxRy(n,:) = nufft_adj(conj(ft_sca_R_y(:,n)).*conj(ft_sca_x).*exp(1i*2*pi*(SamplePoints(:,2))).*wc,stx).';
        end
        x_trunc(a+1:2^R-a, 1:a) = LxMy(2:2^R-2*a+1,1:a);
        x_trunc(a+1:2^R-a, 2^R-a+1:2^R) = RxMy(2:2^R-2*a+1,1:a);
        x_trunc(1:a,a+1:2^R-a) = MxLy(1:a,2:2^R-2*a+1);
        x_trunc(2^R-a+1:2^R, a+1:2^R-a) = MxRy(1:a,2:2^R-2*a+1);

        
        y_concat = FWT2_CDJV_noP(real(x_trunc), ceil(log2(2*a))+1,a)+1i*FWT2_CDJV_noP(imag(x_trunc), ceil(log2(2*a))+1,a);

    end


end


% Copyright (c) 2014. Clarice Poon and Milana Gataric
