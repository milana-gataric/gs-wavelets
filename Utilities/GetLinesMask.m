% Inputs
% N : dimension of mask
% L : number of lines
%
% Output
% mask : N by N matrix with entries 0 or 1. The 1 value entries make up a
% star shape domain of L lines through the centre of the matrix.

function mask = GetLinesMask(N, L)

[x,y] = meshgrid(-floor(N/2):floor(N/2)-1,-floor(N/2):floor(N/2)-1);
mask = zeros(N);
for k=1:L
    ang = k*pi/L;
    mask(y==round(x*cot(ang)))=1;
    mask(x==round(y*tan(ang)))=1;
end

end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
