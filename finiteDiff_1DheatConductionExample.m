%% ============== EnKF - Heat Conduction in a 1D Bar ==================  %%
% Created by Joseph N. Squeo [7/23/2019]
%
% Example taken from:
%[1]        S. Gillijns, O. B. Mendoza, J. Chandrasekar, B. L. R. de Moor,
% D. S. Bernstein, and A. Ridley, What is the ensemble Kalman filter and
% how well does it work?, 2006.
%
% ************************  READ ME  ****************************
% This code utilizes a finite difference method to obtain the temperatures
% at each cell edge of the 1D domain, which can be assimilated into the
% EnKF code as measurements
%% INITIALIZATION 
clear
close all
clc

%---------Material Properties (aluminum)---------%
kTherm = 80; %237;            %[W / m*K] thermal conductivity
Cp = 460; %921;               %[J / kg*K] specific heat @ constant pressure
rho = 7870; %2712;              %[kg / m^3] density
rhoCp = rho * Cp;
alpha = kTherm/rhoCp;    % thermal diffusivity [m^2 / s]

%---------Space---------%
N = 100;            % # cells (since points on cell walls, N+1 points)
cellEdges = 1:N+1;
L = 1;              % length of bar
dx = L/N;           % cell spacing
dx2 = dx*dx;

%---------Time---------%
tEnd = 100;
stabilityCriteria = 0.5*dx2/alpha;   %CFD stability criteria for purely transient diffusion problems
%dt = stabilityCriteria / 5;     
dt = 0.25;
t = 0:dt:tEnd;
tSteps = length(t);
%% =================  FINITE DIFFERENCE METHOD  ======================= %%
%---------Boundary & Initial Conditions, Matrix Initialization---------%
T = zeros(length(t),N+1);   % each row is new time step, each column is new cell edge
Ti = 300;
T(:,1) = Ti;        % boundary condition
T(:,end) = Ti;      % boundary condition
T(1,:) = Ti;        % Initial condition

B = zeros(N+1,1);    % souce input mapping to cells 0.33L and 0.67L
% B( 30:40,1 ) = dt/rhoCp ;
% B( 60:70,1 ) = dt/rhoCp ;
B( round(N*0.33),1 ) = dt/rhoCp ;
B( round(N*0.67),1 ) = dt/rhoCp ;
uk = zeros(1,length(t));   % tEnd+1 because t=0 stored in column 1, t=tEnd stored in column 1+tEnd

%---------Solve for T Using Runge-Kutta Time Integration---------%
q = 5e7;      % W/m^3
uk = q * abs(sin(t));
uk(end+1) = q * abs(sin(tEnd));

for k = 1:length(t)       % TIME LOOP  
    
    for i = 2:N     % SPACE LOOP (points are at cell walls, not center)
        
        qSource = B(i) * uk(k+1);
        
        T(k+1,i) = (dt/rhoCp) * ( (kTherm/dx2) * (T(k,i-1) - 2*T(k,i) ...
            + T(k,i+1)) ) + T(k,i) + qSource;
         
    end
    
end
T(end,:) = [];

% PLOT OF TEMPERATURE DISTRIBUTION VS TIME
figure();
imagesc(cellEdges,t,T);
xlabel('Cell Edge ID')
ylabel('Time (s)')
c = colorbar;
c.Label.String = 'Temperature (K)';
str = sprintf('1D Transient Heat Conduction in a Rod: q_{source}= %.0e*sin(t) W/m^3',q);
title(str);
%% =======================  OPENFOAM  ========================= %%
% fontsize = 17;
% nheader = 22;
% str1 = {};
% str1{end+1} = '~/OpenFOAM/jns14008-5.x/1D_heatConductionEnKF';
% cd(str1{1});
% T = zeros(N,tEnd+1);    % rows are each second of time, columns are cell centers
% T(:,1) = 300;       % initialize temperture [K] at t=0, 0 <= x <= L 
% for time = 1:tEnd
%     str2 = sprintf('%g',time);
%     cd(str2);
%     fid = fopen('T');
%     temp = textscan(fid,'%f', 'headerlines',nheader);
%     fclose(fid);
%     T(:,time+1) = temp{1};
%     cd ../;
% end
% 
% % PLOT OF TEMPERATURE DISTRIBUTION VS TIME
% figure();
% imagesc(T');
% xlabel('Cell Number')
% ylabel('Time (s)')
% c = colorbar;
% c.Label.String = 'Temperature (K)';
% str3 = sprintf('1D Transient Heat Conduction in a Rod: q_{source}= 5e+07 * sin(t) W/m^3');
% title(str3);
%% ==========================  EnKF  =============================== %%
%     Heat conduction in a 1D aluminum rod using the EnKF implementation.    
%     Measurements are available every 10 cells along the 100 cell discretized  
%     rod. The heat source is q = 5000*abs(sin(t)) kW/m^3. T(0,t) = T(L,t)
%     = T(x,0) = 300 K
clc

    %N = 100;                % number of cells
    coeff = dt*alpha/dx2;
    A = zeros(N+1);         % Note the mesh is discretized at the cell edges, not cell centers
    A(1,1) = 1;           % mapping so that T(0,t) = 300K
    A(end,end) = 1;       % mapping so that T(L,t) = 300K
    
    for i = 2:N
        
        % Set up the tridiagonal model operator matrix
        A(i,i-1) = coeff;
        A(i,i) = 1 - (2 * coeff);
        A(i,i+1) = coeff;
    
    end
   
    x_ini = 300 * ones(N+1,1);             % initial model condition T(x,0) = 300K
    sigma = 5;                             % sample error used to generate initial ensemble from x_ini
    w = 2;                               % std dev of model noise
    v = 5;                                 % std dev of measurement noise
    q = 500;                                % # ensemble members (simulation output)
    C = eye(N+1);
    y_meas = nan * ones(N+1,tSteps);   % measurements --> each row = new state, each column = new measurement in time (use nan if no measurement available)
    for cell = 0.1:0.1:0.9
        y_meas(round(N*cell),:) = T(:,round(N*cell));
    end
    %y_meas(:,1) = [];
    [row,col] = size(y_meas);
    tPlot = 10;                         % time (s) for the plot
    
    tic
    [x_EnKF, x_tr,MSE] = ensemblekfilter(A,C,x_ini,w,v,q,y_meas,tSteps,B,uk,sigma,tPlot,dt); 
    toc

   
%% PLOTTING MSE
% For accurate MSE plots, remove random noise from the true state x_tr in
% ensemblekfilter.m

% kArray = 1:k;
% i = 1;
% 
%     for q = [10 20 100]
%         [x_EnKF, x_tr,MSE] = ensemblekfilter(A,C,x_ini,w,v,q,y_meas,tEnd,B,uk,sigma,kplot); 
%         MSE_max(:,i) = max( MSE, [], 2 );  % finds max value across each row to get max MSE value for each time iteration
%         i = i+1;
%     end
%     
% figure;
% box on;
% semilogy(kArray,MSE_max(:,1),':');
% hold on;
% semilogy(kArray,MSE_max(:,2),'-');
% semilogy(kArray,MSE_max(:,3));
% % loglog(kArray,MSE_max(:,4));
% % loglog(kArray,MSE_max(:,5));
% title('EnKF Mean Squared Error vs Ensemble Size');
% ylabel('mean squared error');
% xlabel('time index - k');
% legend('q = 10','q = 20', 'q = 100');
% hold off;

