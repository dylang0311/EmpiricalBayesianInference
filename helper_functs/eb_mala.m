function output = eb_mala(x,mmvMean,params,f_post_ln,lHat,epsilon,sigT,solve_space,m)
% Performs the Metropolis-adjusted Langevin algorithm
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
[N,N_M] = size(x);

num_reject = 0;
num_accept = 0;
output.accept_ratio = zeros(N_M,1);
T = sparse_operator(params);
M = diag(m);

if strcmp(solve_space,'signal')
    A = dftmtx(N);
elseif strcmp(solve_space,'freq')
    A = eye(N);
end

if params.TIMER
    tStart = tic;
end

if params.PRIOR == 1
    mu = 0.0005;
    if strcmp(params.signal,'real')
        y_star = @(x) sign(T*sigT(x)).*min(lHat*abs(T*sigT(x))/mu,1);
    elseif strcmp(params.signal,'complex')
        y_star = @(x) sign(T*abs(sigT(x))).*min(lHat*abs(T*abs(sigT(x)))/mu,1);
    end
    grad_R = @(x) lHat * y_star(x);
else
    grad_R = @(x) 2*(T'*T)*sigT(x);
end

if strcmp(params.signal,'real') && strcmp(solve_space,'signal')
    grad_U = @(x) 2*(real(A'*A)*x-real(A'*mmvMean)) + grad_R(x);
else
    grad_U = @(x) 2*(real(A'*A)*x-A'*mmvMean) + grad_R(x);
end

for kk = 2:N_M
    x_cand = x(:,kk-1) + epsilon^2*M*grad_U(x(:,kk-1))/2+epsilon*sqrt(M)*randn(N,1);
    prop_ln = @(x,y) -1/(2*epsilon^2)*norm(x-y-epsilon^2/2*grad_U(y),2)^2;
    ln_ratio = f_post_ln(x_cand) + prop_ln(x(:,kk-1),x_cand) -...
        f_post_ln(x(:,kk-1)) - prop_ln(x_cand,x(:,kk-1));
    accept_alpha = min(0,ln_ratio);
    u = log(rand);
    if u < accept_alpha
        x(:,kk) = x_cand;
        num_accept = num_accept + 1;
    else
        x(:,kk) = x(:,kk-1);
        num_reject = num_reject + 1;
    end
    output.accept_ratio(kk) = num_accept/(kk-1);
end

output.x = x; 

if params.TIMER
    tEnd = toc(tStart)
end
end

