function [phiNext,num_accept] = malaStepPhi(phiLast,phiPost,num_accept,randNum,xAHy,sigInv)
epsilon = 0.00125;
gradLogPost = @(phi) -epsilon^2*abs(xAHy)*sigInv.*sin(phi-angle(xAHy));
phiCand = phiLast + gradLogPost(phiLast) + epsilon*randNum;
phiCand = mod(phiCand + pi,2*pi) - pi;

prop_ln = @(x,y) -1/(2*epsilon^2)*norm(x-y-gradLogPost(y),2)^2;
ln_ratio = phiPost(phiCand) + prop_ln(phiLast,phiCand) -...
    phiPost(phiLast) - prop_ln(phiCand,phiLast);
accept_alpha = min(0,ln_ratio);
u = log(rand);


if u< accept_alpha
    phiNext = phiCand;
    num_accept = num_accept + 1;
else
    phiNext = phiLast;
end
end