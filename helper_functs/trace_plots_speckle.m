function trace_plots_speckle(x_true,x_mcmc,ind,BI,titlename,fig_name)

figure;
plot(abs(x_mcmc));hold on;
xline(BI,'k--','linewidth',2);
plot(abs(x_true(ind)).*ones(size(x_mcmc)),'r-.','linewidth',2);
xlabel('Iteration');
title(titlename);

if nargin == 5
    saveas(gcf,fig_name)
end

end

