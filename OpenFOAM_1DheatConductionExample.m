
%% ============== EnKF - Heat Conduction in a 1D Bar ==================  %%
% Created by Joseph N. Squeo [Semptember 2019]
%
% Example taken from:
%[1]        S. Gillijns, O. B. Mendoza, J. Chandrasekar, B. L. R. de Moor,
% D. S. Bernstein, and A. Ridley, What is the ensemble Kalman filter and
% how well does it work?, 2006.
%
% ************************  READ ME  ****************************
% This code utilizes OpenFOAM CFD software to obtain the temperatures along
% a 1D mesh. The finite volume method is implemented in OpenFOAM with a
% modified laplacianFoam solver that has the time-step varying source term
% to match the example in the paper cited above. The heat source is q''' =
% abs( 5e7 * sin(t) ) [W/m^3] where t is the time. Note the
% absolute value is used for a heat source rather than a sink.
%% INITIALIZATION 
clear
close all
clc
%---------Material Properties (aluminum)---------%
kTherm = 80;             %[W / m*K] thermal conductivity
Cp = 460;                %[J / kg*K] specific heat @ constant pressure
rho = 7870;              %[kg / m^3] density
rhoCp = rho * Cp;
alpha = kTherm / rhoCp;    % thermal diffusivity [m^2 / s]
%---------Space---------%
N = 100;            % # cells (since points on cell walls, N+1 points)
L = 1;              % length of bar
dx = L/N;           % cell spacing
dx2 = dx*dx;
%---------Time---------%
tEnd = 10;
stabilityCriteria = 0.5*dx2/alpha;   %CFD stability criteria for purely transient diffusion problems
%dt = stabilityCriteria / 5;     
dt = 0.25;
t = 0:dt:tEnd;
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

%% =================  FINITE DIFFERENCE METHOD  ======================= %%
%---------Boundary & Initial Conditions, Matrix Initialization---------%
T_FD = zeros(length(t),N+1);   % each row is new time step, each column is new cell edge
Ti = 300;
T_FD(:,1) = Ti;        % boundary condition
T_FD(:,end) = Ti;      % boundary condition
T_FD(1,:) = Ti;        % Initial condition

B = zeros(N+1,1);    % souce % solverRuns = 1;

%input mapping to cells 0.33L and 0.67L
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
        
        T_FD(k+1,i) = (dt/rhoCp) * ( (kTherm/dx2) * (T_FD(k,i-1) - 2*T_FD(k,i) ...
            + T_FD(k,i+1)) ) + T_FD(k,i) + qSource;
         
    end
    
end
T_FD(end-1:end,:) = [];
T_FD = T_FD(:,1:end-1);

%% ==========================  EnKF  =============================== %%
%     Heat conduction in a 1D aluminum rod using the EnKF implementation.    
%     Measurements are available every 10 cells along the 100 cell discretized  
%     rod. The heat source is q = 5000*abs(sin(t)) kW/m^3. T(0,t) = T(L,t)
%     = T(x,0) = 300 K
clc


x0 = 300;           % initalize ensemble
sigma = 5;          % sample error std deviation to initialize the ensemble
w = 50;            % std dev of model noise
v = 5;              % std dev of measurement noise
q = 30;              % # ensemble members (simulation runs per time step)
%C = eye(N);         % directly maps current state to the forecasted measurements
C = zeros(N);      %%%%%%%%%%%%%%%%%%%%%%%% NEW C MAPPING
solverRuns = 1;
% y_meas = nan * ones(N,tEnd/(dt*solverRuns));   % measurements --> each row = new state, each column = time (Use nan if no measurement available. Num_iterations + 2 because k=1 treated as t=0 & measurement for future state)
y_meas = zeros(N,tEnd/(dt*solverRuns));

% for cell = 0.1:0.1:0.9
%     y_meas(round(N*cell),1:end) = T_FD(solverRuns:solverRuns:end,round(N*cell)); % each column of y_meas is new time, t
%     k = k+1;
% end

for cell = 0.1:0.05:0.9
    C(round(N*cell),round(N*cell)) = 1;  %%%%%%%%%%%%%%%%%%%%%%%% NEW C MAPPING
    y_meas(round(N*cell),1:end) = T_FD(solverRuns:solverRuns:end,round(N*cell)); % each column of y_meas is new time, t
    k = k+1;
end

tPlot = 10;                         % plot at k=10
    if tPlot > tEnd
        error('tPlot (plot time) cannot exceed tEnd (end time)!')
    end
    
% Specify paths to openFOAM case directory and Matlab scripts
solverName = 'myLaplacianFoam';
caseFolder_OF = '~/OpenFOAM/jns14008-5.x/1D_heatConductionEnKF';
scriptsFolder = '~/dataAssimilation';
cd(scriptsFolder);

tic
% [x_EnKF, x_tr] = OF_ensemblekfilter(x0,N,C,w,v,q,y_meas,sigma,tPlot,caseFolder_OF,solverName,tEnd,dt,scriptsFolder,solverRuns);
[x_EnKF, x_tr] = OF_ensemblekfilter2(x0,N,C,w,v,q,y_meas,sigma,tPlot,caseFolder_OF,solverName,tEnd,dt,scriptsFolder,solverRuns);  %%%%%%%%%%%%%%%%%%%%%%%% NEW C MAPPING
toc
x_EnKF1 = x_EnKF;
%%
% solverRuns = 1;
% q = 5;              % # ensemble members (simulation runs per time step)
% y_meas = nan * ones(N,tEnd/(dt*solverRuns));   % measurements --> each row = new state, each column = time (Use nan if no measurement available. Num_iterations + 2 because k=1 treated as t=0 & measurement for future state)
% 
% for cell = 0.1:0.05:0.9
%     y_meas(round(N*cell),1:end) = T_FD(solverRuns:solverRuns:end,round(N*cell)); % each column of y_meas is new time, t
%     k = k+1;
% end
% tPlot = 300;                         % plot at k=10
%     if tPlot > tEnd
%         error('tPlot (plot time) cannot exceed tEnd (end time)!')
%     end
%     
% % Specify paths to openFOAM case directory and Matlab scripts
% solverName = 'myLaplacianFoam';
% caseFolder_OF = '~/OpenFOAM/jns14008-5.x/1D_heatConductionEnKF';
% scriptsFolder = '~/dataAssimilation';
% cd(scriptsFolder);
% 
% tic
% [x_EnKF, x_tr] = OF_ensemblekfilter(x0,N,C,w,v,q,y_meas,sigma,tPlot,caseFolder_OF,solverName,tEnd,dt,scriptsFolder,solverRuns); 
% toc
% x_EnKF2 = x_EnKF;
% %%
% solverRuns = 50;
% q = 5;              % # ensemble members (simulation runs per time step)
% y_meas = nan * ones(N,tEnd/(dt*solverRuns));   % measurements --> each row = new state, each column = time (Use nan if no measurement available. Num_iterations + 2 because k=1 treated as t=0 & measurement for future state)
% 
% for cell = 0.1:0.05:0.9
%     y_meas(round(N*cell),1:end) = T_FD(solverRuns:solverRuns:end,round(N*cell)); % each column of y_meas is new time, t
%     k = k+1;
% end
% tPlot = 300;                         % plot at k=10
%     if tPlot > tEnd
%         error('tPlot (plot time) cannot exceed tEnd (end time)!')
%     end
%     
% % Specify paths to openFOAM case directory and Matlab scripts
% solverName = 'myLaplacianFoam';
% caseFolder_OF = '~/OpenFOAM/jns14008-5.x/1D_heatConductionEnKF';
% scriptsFolder = '~/dataAssimilation';
% cd(scriptsFolder);
% 
% tic
% [x_EnKF, x_tr] = OF_ensemblekfilter(x0,N,C,w,v,q,y_meas,sigma,tPlot,caseFolder_OF,solverName,tEnd,dt,scriptsFolder,solverRuns); 
% toc
% x_EnKF3 = x_EnKF;
% %%
% solverRuns = 100;
% q = 5;              % # ensemble members (simulation runs per time step)
% y_meas = nan * ones(N,tEnd/(dt*solverRuns));   % measurements --> each row = new state, each column = time (Use nan if no measurement available. Num_iterations + 2 because k=1 treated as t=0 & measurement for future state)
% 
% for cell = 0.1:0.05:0.9
%     y_meas(round(N*cell),1:end) = T_FD(solverRuns:solverRuns:end,round(N*cell)); % each column of y_meas is new time, t
%     k = k+1;
% end
% tPlot = 300;                         % plot at k=10
%     if tPlot > tEnd
%         error('tPlot (plot time) cannot exceed tEnd (end time)!')
%     end
%     
% % Specify paths to openFOAM case directory and Matlab scripts
% solverName = 'myLaplacianFoam';
% caseFolder_OF = '~/OpenFOAM/jns14008-5.x/1D_heatConductionEnKF';
% scriptsFolder = '~/dataAssimilation';
% cd(scriptsFolder);
% 
% tic
% [x_EnKF, x_tr] = OF_ensemblekfilter(x0,N,C,w,v,q,y_meas,sigma,tPlot,caseFolder_OF,solverName,tEnd,dt,scriptsFolder,solverRuns); 
% toc
% x_EnKF4 = x_EnKF;
% 
