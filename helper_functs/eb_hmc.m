function output = eb_hmc(x,mmvMean,params,f_post_ln,lHat,epsilon,sigT,solve_space)
% Performs the Hamiltonian Monte Carlo algorithm
%
% Inputs:
%   - x = 2D array of size N1*N2 x N_M
%       - The first column of x is the starting position of the chain
%       - The float values of all columns after the first do not matter
%   - mmvMean = mean of the multiple measurement vectors
%   - params = struct of user-defined parameters
%   - f_post_ln = log posterior density
%   - lHat = point estimate for regularization parameter
%   - epsilon = step size of leapfrog method in HMC
%   - sigT = transform from solve space to signal space
%   - solve_space = space in which we are building the chain
%
% Output:
%   - output = struct with two fields
%       - accept_ratio = acceptance ratio at each step of the chain
%       - x = full chain of size N1*N2 x N_M
%       - hmcEpsilon = final step size of leapfrog method in HMC

[N,N_M] = size(x);

num_reject = 0;
num_accept = 0;
output.accept_ratio = zeros(N_M,1);
leap_steps = 100;
T = sparse_operator(params);

if strcmp(solve_space,'signal')
    A = dftmtx(N);
elseif strcmp(solve_space,'freq')
    A = eye(N);
end

U = @(x) -f_post_ln(x);
if params.PRIOR == 1
    mu = 0.0005;
    if strcmp(params.signal,'real')
        y_star = @(x) sign(T*sigT(x)).*min(lHat*abs(T*sigT(x))/mu,1);
    elseif strcmp(params.signal,'complex')
        y_star = @(x) sign(T*abs(sigT(x))).*min(lHat*abs(T*abs(sigT(x)))/mu,1);
    end
    grad_R = @(x) lHat * y_star(x);
elseif params.PRIOR == 2
    grad_R = @(x) 2*(T'*T)*sigT(x);
end
if strcmp(params.signal,'real') && strcmp(solve_space,'signal')
    grad_U = @(x) 2*(real(A'*A)*x-real(A'*mmvMean)) + grad_R(x);
else
    grad_U = @(x) 2*(real(A'*A)*x-A'*mmvMean) + grad_R(x);
end

for kk = 2:N_M
    [x(:,kk),accept] = hmc_step(U,grad_U,epsilon,leap_steps,x(:,kk-1));
    if accept
        num_accept = num_accept + 1;
    else
        num_reject = num_reject + 1;
    end
    output.accept_ratio(kk) = num_accept/(kk-1);
end

output.x = x; 

end