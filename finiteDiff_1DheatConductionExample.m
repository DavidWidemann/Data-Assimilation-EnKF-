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
L = 1;              % length of bar
dx = L/N;           % cell spacing
dx2 = dx*dx;

%---------Time---------%
tEnd = 1000;
stabilityCriteria = 0.5*dx2/alpha;   %CFD stability criteria for purely transient diffusion problems
%dt = stabilityCriteria / 5;     
dt = 0.25;

%% =================  FINITE DIFFERENCE METHOD  ======================= %%
%---------Boundary & Initial Conditions, Matrix Initialization---------%
T = zeros(tEnd,N+1);   % each row is new time step, each column is new cell edge
Ti = 300;
T(:,1) = Ti;        % boundary condition
T(:,end) = Ti;      % boundary condition
T(1,:) = Ti;        % Initial condition

B = zeros(N+1,1);    % souce input mapping to cells 0.33L and 0.67L
% B( 30:40,1 ) = dt/rhoCp ;
% B( 60:70,1 ) = dt/rhoCp ;
B( round(N*0.33),1 ) = dt/rhoCp ;
B( round(N*0.67),1 ) = dt/rhoCp ;
uk = zeros(1,tEnd);
%---------Solve for T Using Runge-Kutta Time Integration---------%

for k = 1:tEnd       % TIME LOOP
    
    q = 50000;      % W/m^3
    uk(k) = q*10^3 * abs(sin(k));   
    
    for i = 2:N     % SPACE LOOP (points are at cell walls, not center)
        
        qSource = B(i) * uk(k);
        %qSource = 5000e3;    % [W / m^3]
        
        T(k+1,i) = (dt/rhoCp) * ( (kTherm/dx2) * (T(k,i-1) - 2*T(k,i) ...
            + T(k,i+1)) ) + T(k,i) + qSource;
         
    end
    
end

T(k+1,1) = Ti;
T(k+1,end) = Ti;

% PLOT OF TEMPERATURE DISTRIBUTION VS TIME
figure();
imagesc(T);
xlabel('Cell Number')
ylabel('Time (s)')
c = colorbar;
c.Label.String = 'Temperature (K)';
str = sprintf('1D Transient Heat Conduction in a Rod: q_{source}= %d*sin(t) kW/m^3',q);
title(str);

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
    y_meas = nan * ones(N+1,tEnd+1);   % measurements --> each row = new state, each column = new measurement in time (use nan if no measurement available)
    for cell = 0.1:0.1:0.9
        y_meas(round(N*cell),:) = T(:,round(N*cell));
    end
    y_meas(:,1) = [];
    [row,col] = size(y_meas);
    kplot = 10;                         % plot at k=10
    
    [x_EnKF, x_tr,MSE] = ensemblekfilter(A,C,x_ini,w,v,q,y_meas,tEnd,B,uk,sigma,kplot); 
    %[x_EnKF, x_tr,MSE] = HeatConduction1D_ensemblekfilterV3(A,C,x_ini,w,v,q,y_meas,tEnd,B,uk,sigma); 

   
%% PLOTTING MSE
% For accurate MSE plots, remove random noise from the true state x_tr in
% ensemblekfilter.m

kArray = 1:k;
i = 1;

    for q = [10 20 100]
        [x_EnKF, x_tr,MSE] = ensemblekfilter(A,C,x_ini,w,v,q,y_meas,tEnd,B,uk,sigma,kplot); 
        MSE_max(:,i) = max( MSE, [], 2 );  % finds max value across each row to get max MSE value for each time iteration
        i = i+1;
    end
    
figure;
box on;
semilogy(kArray,MSE_max(:,1),':');
hold on;
semilogy(kArray,MSE_max(:,2),'-');
semilogy(kArray,MSE_max(:,3));
% loglog(kArray,MSE_max(:,4));
% loglog(kArray,MSE_max(:,5));
title('EnKF Mean Squared Error vs Ensemble Size');
ylabel('mean squared error');
xlabel('time index - k');
legend('q = 10','q = 20', 'q = 100');
hold off;

