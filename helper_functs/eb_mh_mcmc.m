function output = eb_mh_mcmc(x,params,prop,f_post_ln)
% Performs the Metropolis Hastings MCMC algorithm
%
% Inputs:
%   - x = 2D array of size N1*N2 x N_M
%       - The first column of x is the starting position of the chain
%       - The float values of all columns after the first do not matter
%   - N_M = length of the chain
%   - prop = function handle for proposal density
%   - f_post_ln = log posterior density
%
% Output:
%   - output = struct with two fields
%       - accept_ratio = acceptance ratio at each step of the chain
%       - x = full chain of size N1*N2 x N_M

num_reject = 0;
num_accept = 0;
output.accept_ratio = zeros(params.N_M,1);

if params.TIMER
    tStart = tic;
end

for kk = 2:params.N_M
    x_cand = prop(x(:,kk-1));
    ratio = (f_post_ln(x_cand)-f_post_ln(x(:,kk-1)));
    accept_alpha = min(1,exp(ratio));

    u = rand;

    if u< accept_alpha
        x(:,kk) = x_cand;
        num_accept = num_accept + 1;

    else
        x(:,kk) = x(:,kk-1);
        num_reject = num_reject + 1;
    end

    output.accept_ratio(kk) = num_accept/kk;
end
output.x = x; 

if params.TIMER
    tEnd = toc(tStart)
end
end

