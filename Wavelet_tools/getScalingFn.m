%L : compute scaling function values on k/2^L for k integer
%wname: Daubechies wavelet
function [supp,M] =  getScalingFn(wname, L)
[filter,~,~] = MakeCDJVFilter('LowPass',wname);

% filter = MakeONFilter('Daubechies',2*wname);

N = length(filter);
[X,Y] = meshgrid(1:N,1:N);
X=X(:);
Y = Y(:);

ind = 2*X-Y<1 | 2*X-Y>N;
val_ind = 2*X-Y;
val_ind(ind) = N+1;
filter = [filter,0];
S = sparse(X, Y, filter(val_ind),N,N)*sqrt(2);

val_ind2 = 2*X-Y+1;
ind2 = 2*X-Y+1<1 | 2*X-Y+1>N;
val_ind2(ind2) = N+1;
S2 = sparse(X, Y, filter(val_ind2),N,N)*sqrt(2);

N = length(filter);
[V,E] = eigs(S);
E = diag(E);
[~,ind] = min(abs(E-1));
M = zeros((N-1)*2^L,1);


P = @(x) (x:1:x+N-2)*2^L+1;

vec1 = V(:,ind); %eigenvector corresponding to eigenvalue 1.

M(P(0)) = vec1/sum(vec1); %normalize the eigenvector
V2 = S2*(vec1/sum(vec1));
M(P(1/2)) = V2;

for j=2:L  
    for k=1:2:2^(j-1)-1
        Phi_temp = M(P(k/2^(j-1)));
        
        M(P(k/2^j)) = S*Phi_temp;
        M(P(k/2^j+1/2)) =S2*Phi_temp;

    end
end
supp = (0:(N-1)*2^L-1)/2^L;


end


% Copyright (c) 2014. Clarice Poon and Milana Gataric
