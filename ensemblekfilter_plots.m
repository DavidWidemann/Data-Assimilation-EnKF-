function ensemblekfilter_plots(kplot,q,w,v,y_meas,x_EnKF,x_tr,MSE)

    figure;
    hold on;
    box on;
    plot(y_meas(:,kplot),'kx','linewidth',2);
    plot(x_EnKF(:,kplot),'r:','linewidth',2);
    plot(x_tr(:,kplot),'b:','linewidth',2);
    legend('Observations','EnKF output', 'Truth State w/Noise'); % true value of state by plugging in initial condition to model equations + noise and iterating using previous x_tr
    legend boxoff
    str1 = sprintf('q = %d, k = %d, w = %.1e, v = %.1e', q,kplot,w,v);
    title(str1);
    str2 = sprintf('x_{%d}', kplot);
    ylabel(str2);
    xlabel('Cell Number');
    hold off;

    figure;
    hold on;
    box on;
    [k,n] = size(MSE);
    MSE_max = max( MSE, [], 2 );  % finds max value across each row to get max MSE value for each time iteration
    semilogy(1:k,MSE_max);
    title('EnKF Mean Squared Error vs Ensemble Size');
    ylabel('mean squared error');
    xlabel('time index - k');
    str3 = sprintf('q = %d' , q);
    legend(str3);
    hold off;
    

end