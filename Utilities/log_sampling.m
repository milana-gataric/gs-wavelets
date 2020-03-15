% Generates sampling points of a one dimensional log sampling scheme

% Input: 
%           delta:    delta density of the sampling scheme 
%           K:        a number representing the bandwith, i.e sampling
%                     points are taken from the region [-K,K]
% Output:
%           omega:    one-dimensional vector of length M containing the
%                     sampling points
%           M:        number of sampling points

function [omega,M]=log_sampling(delta,K)

s=0.01;
nu=-log10(delta/2)+s;

K=K+(delta-3*s)/2;

Np=ceil(-(log10(K)+nu)/log10(1-delta/K));  

omega=10.^(-nu+(0:Np)./Np*(log10(K)+nu))-(delta-3*s)/2;
omega=[-fliplr(omega) 0 omega].';

M=length(omega);

end


% Copyright (c) 2014. Clarice Poon and Milana Gataric
