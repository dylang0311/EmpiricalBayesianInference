function [fNoMask,fMask,sigT] = posterior_densities(mmvMean,solve_space,params,...
    lHat,lHatMask,A,AH,M)
% Inputs:
%   - mmvMean = mean of the multiple measurement vectors
%   - solve_space = space in which we are building the chain
%   - params = struct of user-defined parameters
%   - lHat = point estimate for regularization parameter for no mask prior
%   - lHatMask = point estimate for regularization parameter for mask prior
%   - A = forward model
%   - AH = inverse of forward model
%   - M = mask
%
% Outputs:
%   - fNoMask = log posterior density function handle with no mask on prior
%   - fMask = log posterior density function handle with mask on prior
%   - sigT = transform from solve space to signal space


N1 = params.N1;
N2 = params.N2;
L = sparse_operator(params);

if strcmp(solve_space,'signal')
sigT = @(x) x;
    if params.DIM == 1
        if strcmp(params.signal,'real')
            fNoMask = @(x)  (-norm(mmvMean-A(x),2)^2 - ...
                lHat/params.PRIOR * norm(L*x,params.PRIOR)^params.PRIOR);

            fMask = @(x)  (-norm(mmvMean-A(x),2)^2 - ...
                lHatMask/params.PRIOR * norm(M.*L*x,params.PRIOR)^params.PRIOR);
        elseif strcmp(params.signal,'complex')
            fNoMask = @(x)  (-norm(mmvMean-A(x),2)^2 - ...
                lHat/params.PRIOR * norm(abs(L*abs(reshape(x,N1,N2))),params.PRIOR)^params.PRIOR);

            fMask = @(x)  (-norm(mmvMean-A(x),2)^2 - ...
                lHatMask/params.PRIOR * norm(M.*(abs(L*abs(reshape(x,N1,N2)))),params.PRIOR)^params.PRIOR);
        end
    elseif params.DIM == 2
        fNoMask = @(x)  (-norm(reshape(mmvMean-A(x),[],1),2)^2 - ...
            lHat/params.PRIOR * norm(reshape(abs(L*abs(reshape(x,N1,N2))) +...
            abs(abs(reshape(x,N1,N2))*L'),[],1),params.PRIOR)^params.PRIOR);

        fMask = @(x)  (-norm(reshape(mmvMean-A(x),[],1),2)^2 - ...
            lHatMask/params.PRIOR * norm(reshape(M.*(abs(L*abs(reshape(x,N1,N2))) +...
            abs(abs(reshape(x,N1,N2))*L')),[],1),params.PRIOR)^params.PRIOR);
    end
    elseif strcmp(solve_space,'freq')
    if params.DIM == 1
        if strcmp(params.signal,'real')
            sigT = @(x) real(AH(x));
            fNoMask = @(y)  (-norm(mmvMean-y,2)^2 - ...
                lHat/params.PRIOR * norm(freqTV(y,AH,params),params.PRIOR)^params.PRIOR);

            fMask = @(y)  (-norm(mmvMean-y,2)^2 - ...
                lHatMask/params.PRIOR * norm(M.*freqTV(y,AH,params),params.PRIOR)^params.PRIOR);
        elseif strcmp(params.signal,'complex')
            sigT = @(x) AH(x);
            fNoMask = @(y)  (-norm(mmvMean-y,2)^2 - ...
                lHat/params.PRIOR * norm(freqTV(y,AH,params),params.PRIOR)^params.PRIOR);

            fMask = @(y)  (-norm(mmvMean-y,2)^2 - ...
                lHatMask/params.PRIOR * norm(M.*freqTV(y,AH,params),params.PRIOR)^params.PRIOR);
        end
    elseif params.DIM == 2
        if strcmp(params.signal,'real')
            sigT = @(x) real(AH(x));
        elseif strcmp(params.signal,'complex')
            sigT = @(x) AH(x);
        end
        fNoMask = @(y)  (-norm(reshape(mmvMean-y,[],1),2)^2 - ...
            lHat/params.PRIOR * norm(reshape(freqTV(y,AH,params),[],1),params.PRIOR)^params.PRIOR);

        fMask = @(y)  (-norm(reshape(mmvMean-y,[],1),2)^2 - ...
            lHatMask/params.PRIOR * norm(reshape(M.*freqTV(y,AH,params),[],1),params.PRIOR)^params.PRIOR);
    end
end
end