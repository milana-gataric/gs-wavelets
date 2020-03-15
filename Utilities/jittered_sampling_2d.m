% Generates sampling points of a two dimensional jittered sampling scheme
% whose points have density < D in the \ell_1 norm

% Input: 
%           jitter: a number representing jitter
%           D:      a number representing upper bound for the density of
%                   the points in the \ell_1 norm
%           K:      a number representing the bandwith, i.e sampling
%                   points are taken from the region [-K,K]^2
% Output:
%           sp:     two-dimensional vector of sampling points containing
%                   x-coordinate in the first column and y-coordinate in the second
%
% Example:  sp = jittered_samples_2d(0.5,0.1,K);

function sp=jittered_sampling_2d(D,jitter,K)

e=D-jitter-0.0001;
nx=ceil(2*K/e);
e=2*K/nx;

%delta=ni+e;

Kd=K-e;
w=-Kd:e:Kd;
[Dx,Dy]=meshgrid(w,w); 
dx=zeros(length(Dx)); dy=zeros(length(Dx)); 
for k=1:length(Dx)
    dx(k,:)=(-jitter+(2*jitter)*rand(1,length(Dx)));
    dy(k,:)= -(jitter-dx(k,:)) + 2*(jitter-dx(k,:)).*rand(1,length(Dx));
end
Dx=Dx+dx.'; Dy=Dy+dy.'; cx=Dx(:); cy=Dy(:);

ax=-Kd:e:Kd; ay=K*ones(1,length(ax));
dx=(-jitter+(2*jitter)*rand(1,length(ax)));
dy=-abs(-(jitter-dx) + 2*(jitter-dx).*rand(1,length(ax))); 
ax=(ax+dx)'; ay=(ay+dy)';

bx=-Kd:e:Kd; by=-K*ones(1,length(ax));
dx=(-jitter+(2*jitter)*rand(1,length(bx)));
dy=abs(-(jitter-dx) + 2*(jitter-dx).*rand(1,length(bx))); 
bx=(bx+dx)'; by=(by+dy)';

ey=-Kd:e:Kd; ex=-K*ones(1,length(ax));
dx=abs(-jitter+(2*jitter)*rand(1,length(ex)));
dy=-(jitter-dx) + 2*(jitter-dx).*rand(1,length(ex)); 
ex=(ex+dx)'; ey=(ey+dy)';

gy=-Kd:e:Kd; gx=K*ones(1,length(ax));
dx=-abs(-jitter+(2*jitter)*rand(1,length(gx)));
dy=-(jitter-dx) + 2*(jitter-dx).*rand(1,length(gx)); 
gx=(gx+dx)'; gy=(gy+dy)';

sx=[cx; ax; bx; ex; gx; K; K; -K; -K];
sy=[cy; ay; by; ey; gy; K; -K; K; -K];

sp=[sx sy];


end
% Copyright (c) 2014. Clarice Poon and Milana Gataric
