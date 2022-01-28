function ar_plots(ar,BI,fig_name)

figure;
plot(ar,'linewidth',1.5);  hold on;
xline(BI,'k--','linewidth',2);
xlabel('Iteration');
ylim([0,1]);

if nargin == 3
    saveas(gcf,fig_name);
end

end