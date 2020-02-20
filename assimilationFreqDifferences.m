%% Plot differences in frequency plot
Diff1 = abs(x_EnKF1(:,end) - x_EnKF2(:,end));
Diff2 = abs(x_EnKF1(:,end) - x_EnKF3(:,end));
Diff3 = abs( x_EnKF1(:,end) - x_EnKF4(:,end));
err1 = ( (x_EnKF1 - x_EnKF2).^2 );
MSE1 = mean(err1,2);
err2 = ( (x_EnKF1(:,50:50:end) - x_EnKF3(:,1:end-1)).^2 );
MSE2 = mean(err2,2);
err3 = ( (x_EnKF1(:,100:100:end) - x_EnKF4(:,1:end-1)).^2 );
MSE3 = mean(err3,2);
RMSE1 = sqrt(MSE1);
RMSE2 = sqrt(MSE2);
RMSE3 = sqrt(MSE3);

% Plot absolute error differences at latest time
figure(); 
hold on; 
plot(Diff1,'linewidth',2)
plot(Diff2,':','linewidth',2);
plot(Diff3,'--','linewidth',2);
ylabel('Absolute Error');
xlabel('Cell Number');
title('Frequency Error for 1D Heat Conduction Problem');
legend('F=1','F=50','F=100'); legend box off;
box on;
hold off;

m1 = max(Diff1)
m2 = max(Diff2);
m3 = max(Diff3);

% Plot MSE for all times
figure(); 
hold on; 
plot(RMSE1,'linewidth',2)
plot(RMSE2,':','linewidth',2);
plot(RMSE3,'--','linewidth',2);
ylabel('Mean Square Error');
xlabel('Cell Number');
title('Frequency Root-Mean Squared Error');
legend('F=1','F=50','F=100'); legend box off;
box on;
hold off;

maxRMSE1 = max(RMSE1);
maxRMSE2 = max(RMSE2);
maxRMSE3 = max(RMSE3);