function acf_plots(x_mcmc,lag_num,fig_name)

figure;
acf(x_mcmc',lag_num);
xlabel('Lag Length');
if nargin == 3
    saveas(gcf,fig_name)
end

end

