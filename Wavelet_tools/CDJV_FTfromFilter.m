%%computes the fourier transform of Phi at w, \hat Phi(w), where Phi is the LHS scaling
%%function when is_rhs=0, and Phi is the RHS scaling function when
%%is_rhs=1. 
%%iter : the number of iterations
%%LLPEF : filter for the boundary scaling function
%%LPF : filter for the internal scaling function
%%ft0 : fourier transform evaluated at zero, i.e. \hat Phi(0).
%
%%See also, CDJV_FT_Vec.m, CDJV_Setup.m

function FT = CDJV_FTfromFilter(LLPEF, LPF, ft0, w, iter, is_rhs)

N = (length(LLPEF(1,:)) +1) /3; 
ft_phi = FTfromFilter(LPF,w/2,50, 1);

if is_rhs
    m=-N-1:-1:-3*N+1;
    m=m-(N-1);
    ft_phi_vec = ft_phi*exp(-1i*w*m'/2);
else
    m=N:3*N-2; 
    m=m-(N-1);
    ft_phi_vec = ft_phi*exp(-1i*w*m'/2);
end
if iter>0
    ft_interm = CDJV_FTfromFilter(LLPEF, LPF, ft0, w/2,iter-1, is_rhs);
else
    ft_interm = ft0;
end
ft_phiE = [ft_interm; ft_phi_vec];

FT = (1/sqrt(2))*LLPEF * ft_phiE;

% if is_rhs
%     FT = fliplr(FT);
% end

end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
