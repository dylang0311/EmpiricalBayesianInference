function regParams = regularizationParameter(mmv,params,xJumpVec,yJumpVec,M,AH)
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
%   - regParams.lHat = point estimate for regularization parameter for no mask prior
%   - regParams.lHatMask = point estimate for regularization parameter for mask prior


N1 = params.N1;
N2 = params.N2;
L = sparse_operator(params);
mmvMean = 1/params.nMMV*sum(mmv,2); % mean of all the MMVs

data_fidelity_normed = zeros(1,params.nMMV);
reg_normed_laplace = zeros(1,params.nMMV);
reg_normed_SISP = zeros(1,params.nMMV);
for jj = 1:params.nMMV
    data_fidelity_normed(jj) = norm(mmv(:,jj)-mmvMean,2)^2;
%     if strcmp(params.signal,'real')
%         reg_normed_laplace(jj) = norm(xJumpVec(:,jj) + yJumpVec(:,jj),2)^2;
%         reg_normed_SISP(jj) = norm(reshape(M.*reshape(xJumpVec(:,jj) + yJumpVec(:,jj),N1,N2),[],1),2)^2;
%     elseif strcmp(params.signal,'complex')
%         if params.DIM == 1
%             reg_normed_laplace(jj) = norm(L*abs(AH(mmv(:,jj))),2)^2;
%             reg_normed_SISP(jj) = norm(M.*(L*abs(AH(mmv(:,jj)))),2)^2;
%         elseif params.DIM == 2
%             reg_normed_laplace(jj) = norm(L*abs(AH(mmv(:,jj))),2)^2;
%             reg_normed_SISP(jj) = norm(M(:).*L*abs(AH(mmv(:,jj))),2)^2;
%         end
%     end
end

if strcmp(params.signal,'real')
    [xJumpVec,yJumpVec,~] = cFactor(mmvMean,AH,params);
    reg_normed_laplace = norm(xJumpVec + yJumpVec,2)^2;
    reg_normed_SISP = norm(reshape(M.*reshape(xJumpVec + yJumpVec,N1,N2),[],1),2)^2;
elseif strcmp(params.signal,'complex')
    [xJumpVec,yJumpVec,~] = cFactor(mmvMean,AH,params);
    reg_normed_laplace = norm(xJumpVec + yJumpVec,2)^2;
    reg_normed_SISP = norm(reshape(M.*reshape(xJumpVec + yJumpVec,N1,N2),[],1),2)^2;
    if params.DIM == 1
%         reg_normed_laplace = norm(L*abs(AH(mmvMean)),2)^2;
%         reg_normed_SISP = norm(M.*(L*abs(AH(mmvMean))),2)^2;
    elseif params.DIM == 2
%         reg_normed_laplace = norm(L*abs(AH(mmvMean)),2)^2;
%         reg_normed_SISP = norm(M(:).*L*abs(AH(mmvMean)),2)^2;
    end
end

regParams.sigmaSq = 1/(N1*N2)*mean(data_fidelity_normed)/params.nMMV;
regParams.etaSqLaplace = 1/(N1*N2)*mean(reg_normed_laplace);
regParams.etaSqSISP = 1/(nnz(M))*mean(reg_normed_SISP);
if params.PRIOR == 1
    regParams.lHat = 2^(3/2)*regParams.sigmaSq/sqrt(regParams.etaSqLaplace);
    regParams.lHatMask = 2^(3/2)*regParams.sigmaSq/sqrt(regParams.etaSqSISP);
elseif params.PRIOR == 2
    regParams.lHat = regParams.sigmaSq/regParams.etaSqLaplace;
    regParams.lHatMask = regParams.sigmaSq/regParams.etaSqSISP;
end

end