function [lHat,lHatMask] = regularizationParameter(mmv,params,xJumpVec,yJumpVec,M,AH)
% Approximation of regularization parameter lambda done according to
% "Effective new methods for automated parameter selection in
% regularized inverse problems" by Sanders, Platte, and Skeel (2020)
%
% Inputs:
%   - mmv = multiple measurement vectors
%   - params = struct of user-defined parameters
%   - xJumpVec = jump function of first dimension
%   - yJumpVec = jump function of second dimension (zeros if signal 1D)
%   - M = mask
%   - AH = inverse of forward model
%
% Outputs:
%   - lHat = point estimate for regularization parameter for no mask prior
%   - lHatMask = point estimate for regularization parameter for mask prior


N1 = params.N1;
N2 = params.N2;
L = sparse_operator(params);
mmvMean = 1/params.nMMV*sum(mmv,2); % mean of all the MMVs

data_fidelity_normed = zeros(1,params.nMMV);
reg_normed_laplace = zeros(1,params.nMMV);
reg_normed_SISP = zeros(1,params.nMMV);
for jj = 1:params.nMMV
    data_fidelity_normed(jj) = norm(mmv(:,jj)-mmvMean,2)^2;
    if strcmp(params.signal,'real')
        reg_normed_laplace(jj) = norm(xJumpVec(:,jj) + yJumpVec(:,jj),2)^2;
        reg_normed_SISP(jj) = norm(reshape(M.*reshape(xJumpVec(:,jj) + yJumpVec(:,jj),N1,N2),[],1),2)^2;
    elseif strcmp(params.signal,'complex')
        if params.DIM == 1
            reg_normed_laplace(jj) = norm(L*abs(AH(mmv(:,jj))),2)^2;
            reg_normed_SISP(jj) = norm(M.*(L*abs(AH(mmv(:,jj)))),2)^2;
        elseif params.DIM == 2
            reg_normed_laplace(jj) = norm(reshape(abs(L*reshape(abs(AH(mmv(:,jj))),N1,N2))+...
                abs(reshape(abs(AH(mmv(:,jj))),N1,N2)*L.'),[],1),2)^2;
            reg_normed_SISP(jj) = norm(reshape(M.*(abs(L*reshape(abs(AH(mmv(:,jj))),N1,N2))+...
                abs(reshape(abs(AH(mmv(:,jj))),N1,N2)*L.')),[],1),2)^2;
        end
    end
end
sigma_sq = 1/(N1*N2)*mean(data_fidelity_normed);
eta_sq_laplace = 1/(N1*N2)*mean(reg_normed_laplace);
eta_sq_SISP = 1/(nnz(M))*mean(reg_normed_SISP);
if params.PRIOR == 1
    lHat = 2^(3/2)*sigma_sq/sqrt(eta_sq_laplace);
    lHatMask = 2^(3/2)*sigma_sq/sqrt(eta_sq_SISP);
elseif params.PRIOR == 2
    lHat = sigma_sq/eta_sq_laplace;
    lHatMask = sigma_sq/eta_sq_SISP;
end

end