%% Computes the GS reconstruction matrix for 2D Haar wavelets and a uniform 
%  sampling scheme and also the associated minimal singular eigenvalue

clear all;

%% Input Parameters (change as needed)

eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate
S = 1; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 4; %maximum scale of wavelet coefficients

L=2^R/eps;

%% Computation of the GS matrix 

W=(eps/2^R)*(-1*S*L:S*L-1)';
ft_sca= haar_phi_ft(-1*W);


U=zeros((S*2^R*2/eps)^2,2^(2*R));
wc= zeros(2^R,2^R); 

index = 1;


for u=1:2^R
    for v=1:u

        wc(u,v)=1;
        y  = Haar_Op2_Wavelet(1, wc, S, eps, R, ft_sca);
    
        U(:, index) = y(:);
        index=index+1;
        wc(u,v)=0;
        
        if u~= v
            wc(v,u)=1;
            y  = Haar_Op2_Wavelet(1, wc, S, eps, R, ft_sca);
    
            U(:, index) = y(:);
            index=index+1;
            wc(v,u)=0;
        end

       
    end
end

min_sing=sqrt(min(eig(U'*U)));
display(min_sing);

%%
% U_t=zeros(2^(2*R),(S*2^R*2/eps)^2);
% wc= zeros(S*2^R*2/eps,S*2^R*2/eps); 
% 
% index = 1;
% 
% 
% for u=1:2*S*2^R/eps
%     for v=1:u
% 
%         wc(u,v)=1;
%         y  = Haar_Op2_Handle(0, wc, S, eps, R, ft_sca);
%     
%         U_t(:, index) = y(:);
%         index=index+1;
%         wc(u,v)=0;
%         
%         if u~= v
%             wc(v,u)=1;
%             y  = Haar_Op2_Handle(0, wc, S, eps, R, ft_sca);
%     
%             U_t(:, index) = y(:);
%             index=index+1;
%             wc(v,u)=0;
%         end
% 
%        
%     end
% end
% 
% sqrt(min(eig(U_t*U_t')))
% 


% Copyright (c) 2014. Clarice Poon and Milana Gataric
