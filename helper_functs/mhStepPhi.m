function [phiNext,num_accept] = mhStepPhi(phiLast,phiPost,num_accept,randNum)
phiProp = @(phi) mod(normrnd(phi+pi,0.05),2*pi)-pi;
phiCand = phiProp(phiLast);
ratio = phiPost(phiCand) - phiPost(phiLast);

accept_alpha = min(0,ratio);
u = log(randNum);
if u< accept_alpha
    phiNext = phiCand;
    num_accept = num_accept + 1;
else
    phiNext = phiLast;
end
end