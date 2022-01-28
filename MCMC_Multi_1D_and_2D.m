clear all
close all

% This code was originally created by Theresa Scarnati and modified by
% Jiahui Zhang to perform emperical Bayesian inference on 1D real signals
% with additive noise using the Metropolis Hastings MCMC algorithm.

% Support for complex signals, 2D signals, Hamiltonian MC,
% multiplicative noise, and solving in frequency space
% added by Dylan Green.


params.noRNG = 894;
rng(params.noRNG);
%% Parameters
% forward model options
params.signal = 'real'; % 'real' or 'complex' (in complex case, phase random)
params.funct = 'sparseEdge'; % 'sparseSig' or 'sparseEdge'
params.noise_model = 'add'; % 'mult' or 'mult_and_add' or 'add'
params.DIM = 1; % 1 or 2
params.nMMV = 20; % number of MMVs
params.LOOKS = 20; % number of looks (Gamma) (only when noise mult)
params.SNR = 10; % SNR of additive noise (Gaussian) (only when noise add)

% solver options
mcmcSamp = 'metHast'; % 'metHast' or 'hmc'
params.PRIOR = 2; %p of the ell_p params.PRIOR
solveSpace = 'signal'; % 'signal' or 'freq' (freq useful when sparseEdge)

% plotting flags
plot_figs = 1; % diagnosis plots? 
acf_flag = 0; % autocorrelation plots?
trace_flag = 1; % trace plots?
ar_flag = 1;    % acceptance ratio (AR) plots?
postdens_flag = 1;    % posterior density at various pixels plots?

params.N1 = 80; % length of signal dimension 1
N1 = params.N1;
if params.DIM == 1
    params.N2 = 1;
elseif params.DIM == 2
    params.N2 = params.N1; % length of signal dimension 2
end
N2 = params.N2;

params.N_M = 10000; % Chain length
params.BI = 5000; % burn in length

params.COUNTTHRESH = 10; % number of attempts to get the AR acceptable
params.MINRAT = 0.2; % minimal acceptable AR for MH MCMC
params.MAXRAT = 0.8; % maximum acceptable AR for MH MCMC

params.HMCMINRAT = 0.6; % minimal acceptable AR for HMC
params.HMCMAXRAT = 0.95; % minimal acceptable AR for HMC


%% Establish Forward Model, Ground Truth, and MMVs
if params.DIM == 1
    A = @(x) fft(x);
    AH = @(x) ifft(x);
elseif params.DIM == 2
    A = @(x) reshape(fft2(reshape(x,N1,N2,[])),N1*N2,[]);
    AH = @(x) reshape(ifft2(reshape(x,N1,N2,[])),N1*N2,[]);
end

real_and_signal = 0;
if strcmp(params.signal,'real') && strcmp(solveSpace,'signal')
    real_and_signal = 1;
end

[mmv,trueSignal,fig_fold,params,ind,subInd] = build_mmv(params,A);
trueSignalMag = abs(reshape(trueSignal,N1,N2));
if strcmp(params.signal,'complex')
    trueSignalAng = angle(trueSignal);
end

save_file = sprintf('%s_%s_l%d_order%d.png',mcmcSamp,solveSpace,params.PRIOR,params.PAORDER);
mkdir(fig_fold) % adds a folder for figures generated from this execution
mmvMean = 1/params.nMMV*sum(mmv,2); % mean of all the MMVs

%% Concentration Factor Method to get Jump Function if Edges Sparse

if strcmp(params.funct,'sparseEdge')
    xJumpVec = zeros(N1*N2,params.nMMV);
    yJumpVec = zeros(N1*N2,params.nMMV);
    for jj = 1:params.nMMV
        [xJumpVec(:,jj),yJumpVec(:,jj),~] = cFactor(mmv(:,jj),AH,params);
    end
end

%% VBJS
% dimension by dimension
yWeights = zeros(N1,N2).';
yMask = .01*ones(N1,N2).';
if strcmp(params.funct,'sparseEdge')
    xMask = VBJS_weights_speck2D(abs(xJumpVec));
    if params.DIM == 2
        yMask = VBJS_weights_speck2D(abs(yJumpVec));
    end
elseif strcmp(params.funct,'sparseSig')
    xMask = VBJS_weights_speck2D(abs(AH(mmv(:,jj))));
end

M = reshape((xMask+yMask)>1,N1,N2); % mask according to VBJS
m = reshape(M,[],1); % vectorize mask


%% plot signal, indices, and mask
if strcmp(params.signal,'real')
    signalToCompare = trueSignal;
elseif strcmp(params.signal,'complex')
    signalToCompare = trueSignalMag;
end
if params.DIM == 1
    figure;subplot(1,2,1);plot(signalToCompare);hold on
    xline(ind,'k--');title('Magnitude of Signal');
    legend('Signal','Indicators','Location','northwest');hold off
    subplot(1,2,2);plot(M);title('Mask');
elseif params.DIM == 2
    figure;subplot(1,2,1);imagesc(reshape(abs(signalToCompare),N1,N2));
    hold on;plot(subInd(:,1),subInd(:,2),'rx','MarkerSize',15,'LineWidth',2);
    colorbar;title('Magnitude of Signal');hold off
    subplot(1,2,2);imagesc(M);title('Mask');
end
set(gcf,'Position',[100 100 900 400]);
saveas(gcf,[fig_fold,'/originalAndMask_',save_file]);


%% lambda approximations
fprintf('Approximating regularization parameter...\n');

[lHat,lHatMask] = regularizationParameter(mmv,params,xJumpVec,yJumpVec,M,AH);


%% log posterior densities we wish to explore

[f_post_nomask_ln,f_post_mask_ln,sigT] = posterior_densities(mmvMean,solveSpace,params,...
    lHat,lHatMask,A,AH,M);

titlewoutmask = 'Additive Posterior Without Mask';
titlemask = 'Additive Posterior With Mask';

%% MCMC Chains for each Posterior
fprintf('MCMC time...\n');

% starting state of chain
if strcmp(solveSpace,'signal')
    if strcmp(params.signal,'real')
        x_tilde = real(AH(mmvMean));
    elseif strcmp(params.signal,'complex')
        x_tilde = AH(mmvMean);
    end
elseif strcmp(solveSpace,'freq')
    x_tilde = mmvMean;
end

output1 = mcmc_chain(mmvMean,x_tilde,mcmcSamp,solveSpace,params,f_post_nomask_ln,sigT,lHat);
output2 = mcmc_chain(mmvMean,x_tilde,mcmcSamp,solveSpace,params,f_post_mask_ln,sigT,lHatMask);


%% Analysis
ind_1 = ind(1);
ind_2 = ind(2);
ind_3 = ind(3);
ind_4 = ind(4);

indsub = [mod(ind(1),N1),1:N2];


% x1
x1_1 = output1.x(ind_1,:);
x1_2 = output1.x(ind_2,:);
x1_3 = output1.x(ind_3,:);
x1_4 = output1.x(ind_4,:);

% x2
x2_1 = output2.x(ind_1,:);
x2_2 = output2.x(ind_2,:);
x2_3 = output2.x(ind_3,:);
x2_4 = output2.x(ind_4,:);

% error
error = zeros(2,1);
error(1) = norm(trueSignal-mean(output1.x(:,params.BI:end),2))./norm(trueSignal);
error(2) = norm(trueSignal-mean(output2.x(:,params.BI:end),2))./norm(trueSignal);

% calculate relevant means
if strcmp(params.signal,'real')
    output1.image = sigT(mean(output1.x(:,params.BI:end),2));
    output2.image = sigT(mean(output2.x(:,params.BI:end),2));
elseif strcmp(params.signal,'complex')
    output1.image = abs(sigT(mean(output1.x(:,params.BI:end),2)));
    output2.image = abs(sigT(mean(output2.x(:,params.BI:end),2)));

    output1.imageAng = angle(sigT(mean(output1.x(:,params.BI:end),2)));
    output2.imageAng = angle(sigT(mean(output2.x(:,params.BI:end),2)));
end

% calculate 95% confidence intervals
output1c = calc_ci(output1,params.BI,params.signal,sigT); 
output2c = calc_ci(output2,params.BI,params.signal,sigT); 

%% Plot and Save Mean and Credibility Interval Figures
if plot_figs
    fprintf('Calculating and Saving all global figures...\n');
    if params.DIM == 1
        if strcmp(params.signal,'real')
            figure;subplot(1,2,1);
            plot(output1c.image,'b');title('MH No Mask Mean and CI');hold on
            plot(output1c.ci(2,:),'r');plot(output1c.ci(1,:),'r');
            plot(trueSignal,'k');hold off

            subplot(1,2,2);
            plot(output2c.image,'b');title('MH Mask Mean and CI');hold on
            plot(output2c.ci(2,:),'r');plot(output2c.ci(1,:),'r');
            plot(trueSignal,'k');hold off
            saveas(gcf,[fig_fold,'/meanAndCI_',save_file]);
        elseif strcmp(params.signal,'complex')
            figure;subplot(1,2,1);
            plot(output1c.image,'b');title('MH No Mask Magnitude Mean and CI');hold on
            plot(abs(output1c.ci(2,:)),'r');plot(abs(output1c.ci(1,:)),'r');
            plot(trueSignalMag,'k');hold off

            subplot(1,2,2);
            plot(output2c.image,'b');title('MH Mask Magnitude Mean and CI');hold on
            plot(abs(output2c.ci(2,:)),'r');plot(abs(output2c.ci(1,:)),'r');
            plot(trueSignalMag,'k');hold off
            saveas(gcf,[fig_fold,'/meanAndCIMag_',save_file]);

            figure;subplot(1,2,1);
            plot(output1c.imageAng,'b');title('MH No Mask Phase Angle Mean and CI');hold on
            plot(output1c.ciAng(2,:),'r');plot(output1c.ciAng(1,:),'r');
            plot(trueSignalAng,'k');hold off

            subplot(1,2,2);
            plot(output2c.imageAng,'b');title('MH Mask Phase Angle Mean and CI');hold on
            plot(output2c.ciAng(2,:),'r');plot(output2c.ciAng(1,:),'r');
            plot(trueSignalAng,'k');hold off
            saveas(gcf,[fig_fold,'/meanAndCIAng_',save_file]);

            figure;
            phaseErrorLaplace = abs(mod(output1c.imageAng - trueSignalAng - pi/2,pi)-pi/2);
            phaseErrorSISP = abs(mod(output2c.imageAng - trueSignalAng - pi/2,pi)-pi/2);
            plot(phaseErrorLaplace);title(['Phase Error (No Mask avg ',...
                num2str(mean(phaseErrorLaplace)),' and Mask avg ',...
                num2str(mean(phaseErrorSISP)),')']);hold on
            plot(phaseErrorSISP);legend('No Mask','Mask');hold off
            saveas(gcf,[fig_fold,'/phaseError_',save_file]);

            figure;
            plot(abs(output1c.image - trueSignalMag));title('Magnitude Error');hold on
            plot(abs(output2c.image - trueSignalMag));legend('No Mask','Mask');hold off
            saveas(gcf,[fig_fold,'/magnitudeError_',save_file]);
        end

    elseif params.DIM == 2
        figure;subplot(1,2,1);
        imagesc(abs(reshape(output1c.image,N1,N2)));
        colormap gray
        colorbar 
        axis xy image 
        xticks([])
        yticks([])
        title('MH No Mask'); 

        subplot(1,2,2);
        imagesc(abs(reshape(output2c.image,N1,N2)));
        colormap gray
        colorbar 
        axis xy image 
        xticks([])
        yticks([])
        title('MH Mask'); 
        saveas(gcf,[fig_fold,'/mean_',save_file]);

        max_color1 = max(abs(output1c.ci(2,:)-output1c.ci(1,:)),[],'all');
        max_color3 = max(abs(output2c.ci(2,:)-output2c.ci(1,:)),[],'all');
        max_color = max([max_color1 max_color3]);

        figure;subplot(1,2,1);
        imagesc(reshape(abs(output1c.ci(2,:)-output1c.ci(1,:)),N1,N2)); 
        colormap jet; 
        axis xy image 
        xticks([])
        yticks([])
        colorbar
        caxis([0 max_color]);
        title('MH Laplace CI') 


        subplot(1,2,2);
        imagesc(reshape(abs(output2c.ci(2,:)-output2c.ci(1,:)),N1,N2)); 
        colormap jet; 
        axis xy image 
        xticks([])
        yticks([])
        colorbar
        caxis([0 max_color]);
        title('MH SISP CI')
        saveas(gcf,[fig_fold,'/CI_',save_file]);

        if strcmp(params.funct,'sparseSig')
            figure;subplot(1,2,1);
            im1 = reshape(output1c.image,N1,N2);
            ci11 = reshape(output1c.ci(1,:),N1,N2);
            ci12 = reshape(output1c.ci(2,:),N1,N2);
            plot(trueSignalMag(mod(ind(1),N1),1:N2),'k--','LineWidth',1.5);hold on;
            plot(ci11(mod(ind(1),N1),1:N2),'b--','LineWidth',1.5);
            plot(ci12(mod(ind(1),N1),1:N2),'b--','LineWidth',1.5);
            plot(abs(im1(mod(ind(1),N1),1:N2)),'r');hold off;
            ylim([0 2]);xlim([1 N1]);title('No Mask');

            subplot(1,2,2);
            im3 = reshape(output2c.image,N1,N2);
            ci31 = reshape(output2c.ci(1,:),N1,N2);
            ci32 = reshape(output2c.ci(2,:),N1,N2);
            plot(trueSignalMag(mod(ind(1),N1),1:N2),'k--','LineWidth',1.5);hold on;
            plot(ci31(mod(ind(1),N1),1:N2),'b--','LineWidth',1.5);
            plot(ci32(mod(ind(1),N1),1:N2),'b--','LineWidth',1.5);
            plot(abs(im3(mod(ind(1),N1),1:N2)),'r');hold off;
            ylim([0 2]);xlim([1 N1]);title('Mask');

            saveas(gcf,[fig_fold,'/meanAndCISlice_',save_file]);
        elseif strcmp(params.funct,'sparseEdge')
            figure;subplot(1,2,1);
            im1 = reshape(output1c.image,N1,N2);
            ci11 = reshape(output1c.ci(1,:),N1,N2);
            ci12 = reshape(output1c.ci(2,:),N1,N2);
            plot(trueSignalMag(subInd(3,1),1:N2),'k--','LineWidth',1.5);hold on;
            plot(ci11(subInd(3,1),1:N2),'b--','LineWidth',1.5);
            plot(ci12(subInd(3,1),1:N2),'b--','LineWidth',1.5);
            plot(abs(im1(subInd(3,1),1:N2)),'r');hold off;
            ylim([0 2]);xlim([1 N1]);title('No Mask');

            subplot(1,2,2);
            im3 = reshape(output2c.image,N1,N2);
            ci31 = reshape(output2c.ci(1,:),N1,N2);
            ci32 = reshape(output2c.ci(2,:),N1,N2);
            plot(trueSignalMag(subInd(3,1),1:N2),'k--','LineWidth',1.5);hold on;
            plot(ci31(subInd(3,1),1:N2),'b--','LineWidth',1.5);
            plot(ci32(subInd(3,1),1:N2),'b--','LineWidth',1.5);
            plot(abs(im3(subInd(3,1),1:N2)),'r');hold off;
            ylim([0 2]);xlim([1 N1]);title('Mask');
            saveas(gcf,[fig_fold,'/meanAndCISlice_',save_file]);
        end
    end
end

%% Autocorrelation
% reference: https://www.mathworks.com/matlabcentral/fileexchange/30540-autocorrelation-function-acf

if acf_flag
    fprintf('Calculating and Saving all ACF plots...\n')
    
    addpath('./helper_functs/acf/');
    
    lag_num = 5000;
    
    % x1
    acf_plots(abs(x1_1),lag_num,[fig_fold,'/acf_x1_ind1_',save_file]);
    acf_plots(abs(x1_2),lag_num,[fig_fold,'/acf_x1_ind2_',save_file]);
    acf_plots(abs(x1_3),lag_num,[fig_fold,'/acf_x1_ind3_',save_file]);
    acf_plots(abs(x1_4),lag_num,[fig_fold,'/acf_x1_ind4_',save_file]);
    
    
    % x2
    acf_plots(abs(x2_1),lag_num,[fig_fold,'/acf_x2_ind1_',save_file]);
    acf_plots(abs(x2_2),lag_num,[fig_fold,'/acf_x2_ind2_',save_file]);
    acf_plots(abs(x2_3),lag_num,[fig_fold,'/acf_x2_ind3_',save_file]);
    acf_plots(abs(x2_4),lag_num,[fig_fold,'/acf_x2_ind4_',save_file]);
    
end

%% Trace plots

if trace_flag
    fprintf('Calculating and Saving all trace plots...\n');
    % x1
    trace_plots_speckle(trueSignal,x1_1,ind_1,params.BI,'ind_1',...
        [fig_fold,'/trace_x1_ind1_',save_file])
    trace_plots_speckle(trueSignal,x1_2,ind_2,params.BI,'ind_2',...
        [fig_fold,'/trace_x1_ind2_',save_file])
    trace_plots_speckle(trueSignal,x1_3,ind_3,params.BI,'ind_3',...
        [fig_fold,'/trace_x1_ind3_',save_file])
    trace_plots_speckle(trueSignal,x1_4,ind_4,params.BI,'ind_4',...
        [fig_fold,'/trace_x1_ind4_',save_file])
    
    
    % x2
    trace_plots_speckle(trueSignal,x2_1,ind_1,params.BI,'ind_1 M',...
        [fig_fold,'/trace_x2_ind1_',save_file])
    trace_plots_speckle(trueSignal,x2_2,ind_2,params.BI,'ind_2 M',...
        [fig_fold,'/trace_x2_ind2_',save_file])
    trace_plots_speckle(trueSignal,x2_3,ind_3,params.BI,'ind_3 M',...
        [fig_fold,'/trace_x2_ind3_',save_file])
    trace_plots_speckle(trueSignal,x2_4,ind_4,params.BI,'ind_4 M',...
        [fig_fold,'/trace_x2_ind4_',save_file])
    
    
end

%% acceptance ratio plots

if ar_flag
    fprintf('Calculating and Saving all AR plots...\n')
    
    ar_plots(output1c.accept_ratio,params.BI,[fig_fold,'/ar_x1_',save_file])
    ar_plots(output2c.accept_ratio,params.BI,[fig_fold,'/ar_x2_',save_file])
end


%% Probability density over specific points
if postdens_flag
    for ii = 1:length(ind)
        if strcmp(params.signal,'real')
            point_spread1 = output1c.x(ind(ii),params.BI+1:params.N_M);
            point_spread2 = output2c.x(ind(ii),params.BI+1:params.N_M);
            figtitle = sprintf('Probability Density Functions at Index %d',ind(ii));
            signalToCompare = trueSignal;
        elseif strcmp(params.signal,'complex')
            point_spread1 = abs(output1c.x(ind(ii),params.BI+1:params.N_M));
            point_spread2 = abs(output2c.x(ind(ii),params.BI+1:params.N_M));
            figtitle = sprintf('Probability Density Functions of Magnitude at Index %d',ind(ii));
            signalToCompare = trueSignalMag;
        end

        subplot(2,2,ii);kde(point_spread1);hold on
        kde(point_spread2);
        xline(signalToCompare(ind(ii)),'--k','LineWidth',2);
        title(figtitle);
        legend('No Mask','Mask','True Value');
        set(gcf,'Position',[100 100 1000 800]);
        saveas(gcf,[fig_fold,'/probDensities_',save_file]);
    end
end
