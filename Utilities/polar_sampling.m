% Generates sampling points of a two dimensional polar (radial) sampling scheme
% whose points have density < D in the \ell_2 norm
%
% Input: 
%           D:      a number representing the upper bound for the density of
%                   the points in the \ell_2 norm
%           K:      a number representing the bandwith, i.e sampling
%                   points are taken from the Euclidean ball of radius K
% Output:
%           sp:     two-dimensional vector of sampling points containing
%                   x-coordinate in the first column and y-coordinate in the second
%           nr:     the number of radial lines along which sampling points are taken 
%
% Example:  [sp, nr] = polar_sampling(K,sqrt(2)/4);

function [sp, nr]=polar_sampling(K, D)

dr = 0.4; %dinstance on a radial line; 

dt=2*(atan(sqrt(D^2-(dr/2)^2)/(K-dr/2))); 
nr=ceil(pi/dt)+1; % number of radial lines such that delta_\ell^2(sp) < D
dtheta=pi/nr;

theta=0:dtheta:pi-dtheta;
r=-K:dr:K;
sx=[]; sy=[];
for j=1:length(theta)
    sx=[sx r*cos(theta(j))]; %#ok<AGROW>
    sy=[sy r*sin(theta(j))]; %#ok<AGROW>
end
Us=unique([sx' sy'],'rows'); sx=Us(:,1)'; sy=Us(:,2)';

sp(:,1)=sx.'; sp(:,2)=sy.';


end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
