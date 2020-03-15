% Returns a mask of length M containing k entries of value 1, and the remaining entries of value 0. 
% By considering the mask to be indexed by -M/2+1:M/2, we divide the indices into Lev levels, 
% where index k is in the n^{th} level are indices such that M*(n-1)/(2*Lev)<=|k|<M*n/(2*Lev).
% Then, if index k belongs to level n,
% P( mask(k)=1) = C exp(-(bn/N)^a), for some appropriate constant C such that this is a probability distribution.
% In addition, we let mask(k)=1 for all k such that |k|<=c.


%Inputs:
% M: size of the sampling mask
% k: number of entries set to one
% a,b: parameters governing the distribution which we draw from, 
%       larger values of a and b will mean a higher concentration on the centre of the mask.
% c: specifies the entries which will be set to one
% Lev: number of levels.

%Output:
% mask: M vector whose entries are either zero or one.

function mask = GetSampleMask1D(M, k, a,b,c,Lev)
mask = zeros(M,1);
total=0;
upper = max(k-c,1);

indx = -M/2+1:M/2;
centre = find(abs(indx)<=c);

N = M/2;
w = exp(-(b*ceil(abs(indx)*Lev/N)/Lev).^a);
w = w(:);


while total<k
    lower = upper;
    upper = 2*upper;
    samples = randsample(M, upper,  true, w);
    mask(samples)=1; mask(centre)=1;
    total = sum(mask);
    mask(samples)=0; mask(centre)=0;
end

while upper-lower>1
    mid = ceil((upper+lower)/2 );
    samples = randsample(M, upper,  true, w);
    mask(samples)=1;  mask(centre)=1;
    total = sum(mask);
    mask(samples)=0; mask(centre)=0;
    if total>k
        upper = mid;
    else
        lower = mid;
    end
end

mask(samples)=1; mask(centre)=1;

% Copyright (c) 2014. Clarice Poon and Milana Gataric
