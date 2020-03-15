% Returns a sampling mask with k entries of value 1, and the remaining entries of value 0. 
% By considering the matrix to be indexed by -M/2+1:M/2 along the horizontal and vertical axis, 
% we divide the indices into Lev levels via concentric circles centred at the index (0,0).
% If index (k,j) belongs to level n, P( mask(k,j)=1) = C exp(-(bn/N)^a),
%  for some appropriate constant C such that this is a probability distribution.
% In addition, we let mask(k,j)=1 for all (k,j) such that |(k,j)|<=c.


%Inputs:
% M: size of the sampling mask
% k: number of entries set to one
% a,b: parameters governing the distribution which we draw from, larger values of a and b will mean a higher concentration on the centre of the mask.
% c: specifies the entries which will be set to one
% Lev: number of levels.

%Output:
% mask: M by M matrix whose entries are either zero or one.




function mask = GetSampleMask2D(M, k, a,b,c,Lev)
mask = zeros(M);
total=0;
upper = max(k-c^2,1);

[xInd,yInd] = meshgrid(-M/2+1:M/2,-M/2+1:M/2);
centre = find(sqrt(xInd.^2+yInd.^2)<=c);
cInd = xInd+1i*yInd;
N = max(abs(cInd(:)));

w = exp(-(b*ceil(abs(cInd)*Lev/N)/Lev).^a);
w = w(:);


while total<k
    lower = upper;
    upper = 2*upper;
    samples = randsample(M^2, upper,  true, w);
    mask(samples)=1; mask(centre)=1;
    total = sum(mask(:));
    mask(samples)=0; mask(centre)=0;
end

while upper-lower>1
    mid = ceil((upper+lower)/2 );
    samples = randsample(M^2, upper,  true, w);
    mask(samples)=1;  mask(centre)=1;
    total = sum(mask(:));
    mask(samples)=0; mask(centre)=0;
    if total>=k
        upper = mid;
    else
        lower = mid;
    end
end

mask(samples)=1;  mask(centre)=1;
end



% Copyright (c) 2014. Clarice Poon and Milana Gataric
