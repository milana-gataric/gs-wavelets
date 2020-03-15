% Generates sampling points of a one dimensional sampling scheme

% Input: 
%           nu:         a number representing jitter
%           epsilon:    a number representing fixed part of the distance
%                       between the points
%           K:          a number representing the bandwith, i.e sampling
%                       points are taken from the region [-K,K]
% Output:
%           omega:      one-dimensional vector of length M containing the
%                       sampling points
%           M:          number of sampling points

function [omega,M]=jittered_sampling_1d(nu,epsilon,K) 

M=ceil((2*K-nu)/epsilon)+1; 
omega=zeros(M,1);
omega(1)=-K;
omega(2:M-1)=-K+(1:M-2)*epsilon+(-nu+(2*nu)*rand(1,M-2));
omega(M)=K;

end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
