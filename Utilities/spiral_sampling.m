% Generates sampling points of a two dimensional spiral sampling scheme
% which can be a single spiral or have more than 1 leaves and  whose points 
% have density < D in the \ell_2 norm
%
% Input: 
%           D:      a number representing the upper bound for the density of
%                   the points in the \ell_2 norm
%           K:      a number representing the bandwith, i.e sampling
%                   points are taken from the Euclidean ball of radius K
%           nl:     an integer representing number of spiral's leaves
% Output:
%           sp:     two-dimensional vector of sampling points containing
%                   x-coordinate in the first column and y-coordinate in the second 
%
% Example:  sp = spiral_sampling(K,sqrt(2)/4,4);


function sp = spiral_sampling(K, D, nl)

% spiral with constant angular velocity
dr=nl*0.5; %dinstance on a radial line
k=K/dr;

if  nl==1


    xfun =  @(t,k) ((k-t/(2*pi)).^2+k*(k-1)-cos(t)*(2*k-1)*(k-t/(2*pi)))/(sin(t)*(2*k-1)*(k-t/(2*pi)));
    dfun = @(t,c,k) (c/2)^2 + (c*k-c/2).^2 * xfun(t,k).^2  - D^2;                  
    fun = @(x) dfun(x,dr,k);  
    t0 = [0.001,pi-0.001];
    dt = fzero(fun,t0);

    nr=ceil((2*pi)/dt)+1;
    dt=2*pi/(nr);

    theta=0:dt:(k*2*pi);
    sx=dr/(2*pi)*theta.*real(exp(1i*theta));
    sy=dr/(2*pi)*theta.*imag(exp(1i*theta));
    
    
elseif  nl>1

    xfun =  @(t,k) ((k-t/(2*pi)).^2+k*(k-1)-cos(t)*(2*k-1)*(k-t/(2*pi)))/(sin(t)*(2*k-1)*(k-t/(2*pi)));
    dfun = @(t,c,k) ((c/nl)/2)^2 + (c*k-c/2).^2 * xfun(t,k).^2  - D^2;                  
    fun = @(x) dfun(x,dr,k);  
    t0 = [0.001,pi-0.001];
    dt = fzero(fun,t0);

    nr=ceil((2*pi)/dt)+1;
    dt=2*pi/(nr);

    sx=[];sy=[];
    theta=0:dt:k*2*pi;
    for m=1:4
        sx=dr/(2*pi)*theta.*real(exp(1i*(theta-2*pi*m/nl)));
        sy=dr/(2*pi)*theta.*imag(exp(1i*(theta-2*pi*m/nl)));
    end

end

sp(:,1)=sx.'; sp(:,2)=sy.';

end


% Copyright (c) 2014. Clarice Poon and Milana Gataric
