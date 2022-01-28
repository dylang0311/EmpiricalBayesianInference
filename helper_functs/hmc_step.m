function [return_q,accept] = hmc_step(U,grad_U,epsilon,L,current_q)
% Adapted from MCMC Handbook Chapter 5.
% Performs the HMC update step

q = current_q;
p = normrnd(0,1,length(q),1);
current_p = p;

% Make a half step for momentum at the beginning
p = p - epsilon * grad_U(q) / 2;

% Alternate full steps for position and momentum
for ii = 1:L
    % Make a full step for the position
    q = q + epsilon * p;
    % Make a full step for the momentum, except at end of trajectory
    if ii ~= L
        p = p - epsilon * grad_U(q);
    end
end
% Make a half step for momentum at the end
p = p - epsilon * grad_U(q) / 2;
% Negate momentum at end of trajectory to make the proposal symmetric
p = -p;

% Evaluate potential and kinetic energies at start and end of
% trajectory

current_U = U(current_q);
current_K = sum(current_p.^2)/2;
proposed_U = U(q);
proposed_K = sum(p.^2)/2;

% Accept or reject the state at tend of trajectory, returning either
% the position at the end of the trajectory or the initial position

u = rand;
test = current_U-proposed_U+current_K-proposed_K;
if u < exp(test)
    return_q = q; % accept
    accept = 1;
else
    return_q = current_q; % reject
    accept = 0;
end
end