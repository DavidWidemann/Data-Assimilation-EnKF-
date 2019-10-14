function MSE_plots(MSE,q) 
%% Plot MSE if passed through to the function
    figure;
    hold on;
    box on;
    [k,n] = size(MSE);
    MSE_max = max( MSE, [], 2 );  % finds max value across each row to get max MSE value for each time iteration
    semilogy(1:k,MSE_max);
    title('EnKF Mean Squared Error vs Ensemble Size');
    ylabel('mean squared error');
    xlabel('time index - k');
    str = sprintf('q = %d' , q);
    legend(str);
    hold off;
end