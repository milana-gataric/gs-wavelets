function ft = FTfromFilter(Lo_R,x,iter, isScalingFn)

if isScalingFn
    p=1:iter;
    ft = prod(hat_h(x./2.^p,Lo_R));
else
    ft = -1*exp(-1i*x/2)*conj(hat_h(x/2+pi,Lo_R))*FTfromFilter(Lo_R,x/2,iter, 1);
end

end

function H = hat_h(w,h)
H=0;
for N= 0:length(h)-1
    H=H+ h(N+1).*exp(-1i*N.*w);
end
H=H./sqrt(2);
end



% Copyright (c) 2014. Clarice Poon and Milana Gataric
