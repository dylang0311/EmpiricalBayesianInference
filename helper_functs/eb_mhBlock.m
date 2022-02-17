function output = eb_mhBlock(x,params,m,regParams,mmvMean,A,AH,etaInv)
% Performs the Metropolis Hastings MCMC algorithm
%
% Inputs:
%   - x = 2D array of size N1*N2 x N_M
%       - The first column of x is the starting position of the chain
%       - The float values of all columns after the first do not matter
%   - N_M = length of the chain
%   - prop = function handle for proposal density
%   - f_post_ln = log posterior density
%   - m = mask
%
% Output:
%   - output = struct with two fields
%       - accept_ratio = acceptance ratio at each step of the chain
%       - x = full chain of size N1*N2 x N_M

num_reject = 0;
num_accept = 0;
output.accept_ratio = zeros(params.N_M,1);

[T] = build_mask(m,params);
N1 = params.N1;
N2 = params.N2;
    

sigInv = 1/regParams.sigmaSq;

if params.TIMER
    tStart = tic;
end

if strcmp(params.signal,'real')
    invMat = inv(etaInv*(T'*T)+sigInv*N1*N2*eye(N1*N2));
    xMean = invMat*(sigInv*N1*N2*real(AH(mmvMean)));
    x = mvnrnd(xMean,invMat,params.N_M).';
    output.x = x; 
    
elseif strcmp(params.signal,'complex')
    mag = abs(x);
    theta = @(phi) exp(1i*(phi));
    thetaConj = @(phi) exp(-1i*(phi));
    AHy = AH(mmvMean);
    
    phi = x;
    phi(:,1) = angle(x(:,1));

    
    invMat = inv(etaInv*(T'*T)+N1*N2*2*sigInv*eye(N1*N2));
    invMatAHy = invMat*2*sigInv*N1*N2*AHy;
    xVarChol = chol(invMat);
    ranNum = randn(N1*N2,params.N_M);
%     ranNum2 = rand(1,params.N_M);
    ranNum2 = randn(N1*N2,params.N_M);
    for kk = 2:params.N_M
        xMean = invMat*real(2*sigInv*thetaConj(phi(:,kk-1))*N1*N2.*AHy);
        currentRanNum = ranNum(:,kk);
        while sum(xMean + xVarChol*currentRanNum <= 0) > 0
            currentRanNum = randn(N1*N2,1);
        end
        mag(:,kk) = xMean + xVarChol*currentRanNum;
        
        xAHy = mag(:,kk).*AHy;
        phiPost = @(phi) sum(sigInv*2*N1*N2*abs(xAHy).*cos(phi-angle(xAHy)));
        [phi(:,kk),num_accept] = malaStepPhi(phi(:,kk-1),phiPost,num_accept,ranNum2(:,kk),xAHy,sigInv);
        output.accept_ratio(kk) = num_accept/kk;   
    end
    
    output.x = mag.*exp(1i*phi); 
    
end


if params.TIMER
    tEnd = toc(tStart)
end
end

