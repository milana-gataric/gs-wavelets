%% Generates 2D CDJV matrix when samples are taken uniformly

clear all; 

%% Input Parameters (change as needed)

a = 2; %number of vanishing moments; it can be in {1,2,3}
eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 4; %maximum scale of wavelet coefficients

L = 2^R/eps;

%% Compute the fourier transforms of the scaling functions

disp('Precomputing Fourier transforms of the scaling functions...')

[ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);

%% Compute the matrix

disp('Computing the GS matrix...')
U=zeros((S*L*2)^2,2^(2*R));
wc= zeros(2^R,2^R); 

index = 1;

for u=1:2^R
    for v=1:u

        wc(u,v)=1;
        y  = CDJV_Op2_Wavelet(1, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
    
        U(:, index) = y(:);
        index=index+1;
        wc(u,v)=0;
        
        if u~= v
            wc(v,u)=1;
            y  = CDJV_Op2_Wavelet(1, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
    
            U(:, index) = y(:);
            index=index+1;
            wc(v,u)=0;
        end

       
    end
end

min_sing = sqrt(min(eig(U'*U)));
display(min_sing)

%% Transpose
% % NB: Ut is a matrix representation of the adjoint operator, but not the
% % adjoint of the matrix U above, because I'm not being consistent in the ordering of
% % the vectors.
% 
% Ut=zeros(2^(2*R),(S*L*2)^2);
% wc= zeros(S*L*2,S*L*2); 
% 
% index = 1;
% 
% 
% for u=1:S*L*2
%     for v=1:u
% 
%         wc(u,v)=1;
%         y  = CDJV_Op2_Wavelet(0, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
%     
%         Ut(:, index) = y(:);
%         index=index+1;
%         wc(u,v)=0;
%         
%         if u~= v
%             wc(v,u)=1;
%             y  = CDJV_Op2_Wavelet(0, wc, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R);
%     
%             Ut(:, index) = y(:);
%             index=index+1;
%             wc(v,u)=0;
%         end
% 
%        
%     end
% end
% 
% sqrt(min(eig(Ut*Ut')))


% Copyright (c) 2014. Clarice Poon and Milana Gataric
