function output = mcmc_chain(mmvMean,x_tilde,mcmc_samp,solve_space,params,logPosterior,sigT,lHat)
% performs the selected MCMC method for sampling from the posterior
%
% Inputs:
%   - mmvMean = mean of the multiple measurement vectors
%   - x_tilde = starting position of chain
%   - mcmc_samp = MCMC method
%   - solve_space = space in which we are building the chain
%   - params = struct of user-defined parameters
%   - logPosterior = log of posterior we wish to sample from
%   - sigT = transform from solve space to signal space
%   - lHat = point estimate for regularization parameter
%
% Output:
%   - output = struct with two or three fields
%       - x = full MCMC chain
%       - accept_ratio = acceptance ratio at each step of the chain
%       if mcmc_samp = 'hmc'
%           - hmcEpsilon = final step size of leapfrog method in HMC


N1 = params.N1;
N2 = params.N2;
output.x = zeros(N1*N2,params.N_M);
output.x(:,1) = x_tilde; % initial condition

real_and_signal = 0;
if strcmp(params.signal,'real') && strcmp(solve_space,'signal')
    real_and_signal = 1;
end

switch mcmc_samp
    case{'metHast'}
    
    % initial proposal distribution variance
    propStd = 0.5;
    
    % initial proposal distribution (Gaussian)
    if real_and_signal
        prop = @(x) normrnd(x,propStd);
    else
        prop = @(x) (normrnd(real(x),propStd) +...
            1i*normrnd(imag(x),propStd));
    end
    output.accept_ratio = 0; % initialize acceptance ratio
    cnt = 0; % initialize count of times MCMC chain generated

    while or(mean(output.accept_ratio) < params.MINRAT, mean(output.accept_ratio) > params.MAXRAT)
        
        % run chain
        output = eb_mh_mcmc(output.x,params.N_M,prop,logPosterior);
        
        % update variances
        if mean(output.accept_ratio) < params.MINRAT
            propStd = propStd/2;
        elseif mean(output.accept_ratio) > params.MAXRAT
            propStd = propStd*2;
        end
        
        % update proposal densities
        if real_and_signal
            prop = @(x) normrnd(x,propStd);
        else
            prop = @(x) (normrnd(real(x),propStd) +...
                1i*normrnd(imag(x),propStd));
        end

        % check if number of iterations exceeds threshold
        cnt = cnt+1;
        if cnt > params.COUNTTHRESH
            fprintf('Unable to reach appropriate acceptance ratio in posterior # 1 \n')
            break
        end
    end
    output.propStd = propStd;
    
    case 'hmc'
    hmcEpsilon = 0.001; % initialize stepsize for leapfrog solver
    output = eb_hmc(output.x,mmvMean,params,logPosterior,lHat,hmcEpsilon,sigT,solve_space);
    cnt = 1;
    while or(mean(output.accept_ratio) < params.HMCMINRAT, mean(output.accept_ratio) > params.HMCMAXRAT)
        if cnt > params.COUNTTHRESH
            fprintf('Unable to reach appropriate acceptance ratio in posterior # 1 \n')
            break
        end
        if mean(output.accept_ratio(params.BI:end)) > params.HMCMINRAT
            hmcEpsilon = hmcEpsilon * 2;
        elseif mean(output.accept_ratio(params.BI:end)) < params.HMCMAXRAT
            hmcEpsilon = hmcEpsilon / 2;
        end
        output = eb_hmc(output.x,mmvMean,params,logPosterior,lHat,hmcEpsilon,sigT,solve_space);
        cnt = cnt+1;
    end
    output.hmcEpsilon = hmcEpsilon;
        
end

output.x = reshape(output.x,N1*N2,params.N_M);
end