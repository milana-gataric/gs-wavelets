% Output: 
%           ft: Fourier transform of Harr scaling function evaluated at 
%               the points specified in the vector x
% Input:
%           x:  A vector of frequences where the Fourier transform is evaluated 
%

function ft =haar_phi_ft(x)

ft = zeros(length(x),1);

for j=1:length(x)
    if(x(j)==0)
        ft(j)=1;
    else
        ft(j) = (exp(-1i*2*pi*x(j))-1)./(-2*1i*pi*x(j));
    end
end

end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
