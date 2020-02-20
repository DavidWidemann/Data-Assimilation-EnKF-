function ensemblekfilter_plots(tPlot,q,w,v,y_meas,x_EnKF,x_tr,dt,solverRuns)
   kPlot = tPlot/(dt*solverRuns);    % Index for specified time to access matrices
    %% Plot observations, EnKF output and true state
    figure(1);
    hold on;
    box on;
    plot(y_meas(:,kPlot),'kx','linewidth',2);
    plot(x_EnKF(:,kPlot),'r-','linewidth',2);
    plot(x_tr(:,kPlot),'b:','linewidth',2);
    legend('Observations','EnKF output', 'True State'); % true value of state by plugging in initial condition to model equations + noise and iterating using previous x_tr
%% Frequency Tests
%     plot(y_meas(:,kPlot),'x','linewidth',2);
%     plot(x_EnKF(:,kPlot),'--','linewidth',2);
%     legend('Meas','Original k=80W/m*K','Meas F=1','EnKF F=1', 'Meas F=50','EnKF F=50','Meas F=100','EnKF F=100');
%     ylim([200 700]);  
%%
    legend boxoff
    str1 = sprintf('q = %g, t = %g s, w = %g, v = %g', q,tPlot,w,v);
    title(str1);
    str2 = sprintf('T(t = %g s)', tPlot);
    ylabel(str2);
    xlabel('Cell Number');
    ylim([200 400]);
    hold off;
    
end
