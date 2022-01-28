function [mmv,trueSignal,fig_fold,params,ind,subInd] = build_mmv(params,A)
% Inputs:
%   - params = struct of user-defined parameters
%   - A = forward model
%
% Outputs:
%   - mmv = multiple measurement vectors
%   - trueSignal = the true underlying signal
%   - fig_fold = string; name of the figure folder to be made
%   - params = struct of user-defined parameters with additional field
%       - PAORDER = order of polynomial annihilation operator needed
%   - ind = indices to be used for analysis
%   - subInd = if dim=2 and sparse in edges, subscripts of ind; else -1

N1 = params.N1;
N2 = params.N2;

a = 0; % lower bound of domain
b = 1; % upper bound of domain
dx = (b-a)./(N1-1); % interval size
spatGrid = a+ dx*(0:N1-1)'; % gridpoints

if strcmp(params.funct,'sparseSig')
    s = 4; % number of nonzero pixels
    g = zeros(N1,N2); 
    ind = randperm(N1*N2,s); 
    g(ind) = 1; 
    
    params.PAORDER = 0; % Polynomial Annihilation (PA) order
elseif strcmp(params.funct,'sparseEdge')
    form = 'sqCirc';
    if params.DIM == 1
        g = zeros(N1,1);
        g(0.1 < spatGrid & spatGrid < 0.25) = 40;
        g(0.32 < spatGrid & spatGrid < 0.35) = 10;
        g(0.6 <= spatGrid & spatGrid <= 0.9) =...
            50*sin(10/3*pi*spatGrid(0.6 <= spatGrid & spatGrid <= 0.9)).^2;
        g = (g+50)./100;
        ind = [find(spatGrid > 0.1, 1); find(spatGrid > 0.32, 1);
                find(spatGrid >= 0.6, 1); find(spatGrid >= 0.75, 1)];
        subInd = -1;
    elseif params.DIM == 2
        [g,ind,subInd] = generate_test_image(form,N1,N2);
    end
    
    params.PAORDER = 1; % PA order
end

if strcmp(params.signal,'real')
    trueSignal = g;
    trueSignal = reshape(trueSignal,N1*N2,1);
elseif strcmp(params.signal,'complex')
    phi = 2*pi*(rand(N1,N2)-0.5);
    Theta = exp(1i*phi);
    trueSignal = Theta.*g;
    
    trueSignal = reshape(trueSignal,N1*N2,1);
end

noiseStandardDeviation = mean(abs(A(trueSignal)),'all')*10^(-(params.SNR/20));

nu = gamrnd(params.LOOKS,1/params.LOOKS,N1*N2,params.nMMV); % multiplicative noise
eta = noiseStandardDeviation/sqrt(2)*(randn(N1*N2,params.nMMV) +...
    1i*randn(N1*N2,params.nMMV)); % additive noise
mmv = zeros(N1*N2,params.nMMV); % initialize MMVs

switch params.noise_model
    case 'add'
        mmv = A(trueSignal) + eta;
        fig_fold = sprintf('./Figures/results_%s_%s_%dD_%dSNR_%dMMV_RNG%d',...
            params.signal,params.funct,params.DIM,params.SNR,params.nMMV,params.noRNG);
    case 'mult'
        mmv = A(nu.*trueSignal);
        fig_fold = sprintf('./Figures/results_%s_%s_%dD_%dLooks_%dMMV_RNG%d',...
            params.signal,params.funct,params.DIM,params.LOOKS,params.nMMV,params.noRNG);
    case 'mult_and_add'
        mmv = A(nu.*trueSignal) + eta;
        fig_fold = sprintf('./Figures/results_%s_%s_%dD_%dLooks_%dSNR_%dMMV_RNG%d',...
            params.signal,params.funct,params.DIM,params.LOOKS,params.SNR,params.nMMV,noRNG);
end

end