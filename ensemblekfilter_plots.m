
function ensemblekfilter_plots(tPlot,q,w,v,y_meas,x_EnKF,x_tr,dt)
   kPlot = tPlot/dt;    % Index for specified time to access matrices
    %% Plot observations, EnKF output and true state
    figure;
    hold on;
    box on;
    plot(y_meas(:,kPlot),'kx','linewidth',2);
    plot(x_EnKF(:,kPlot),'r:','linewidth',2);
    plot(x_tr(:,kPlot),'b:','linewidth',2);
    legend('Observations','EnKF output', 'True State'); % true value of state by plugging in initial condition to model equations + noise and iterating using previous x_tr
    legend boxoff
    str1 = sprintf('q = %g, t = %g s, w = %g, v = %g', q,tPlot,w,v);
    title(str1);
    str2 = sprintf('T(t = %g s)', tPlot);
    ylabel(str2);
    xlabel('Cell Number');
    hold off;
    
end
