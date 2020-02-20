
%   ===================================================================   %
%   Created by Joseph N. Squeo, University of Connecticut
%   September 2019
%   ===================================================================   %
%	
%   The ensemblekfilter function is an ensemble kalman filter used to
%	predict a future state in time of a system. I have written a document
%	that corresponds to the numbered steps in this file. The equations are
%	explained and the purpose of each step to provide a better
%	understanding of the code and how it works. 
%
%   There are two states: truth state (x_tr) and estimate state (x_est).
%   The ensemble kalman filter initially generates an ensemble of states by
%   adding random sample error to the initial state. The initial state and
%   standard deviation of the random sample error are function inputs.
%	x_(k+1) = f(x,u) + w, where u some input, w the Gaussian distributed
%	process noise, and f is a nonlinear function. The measurement is
%	y_(k+1) = h(x) + v where h is a nonlinear function and v Gaussian
%	distributed measurement noise.
%   
%   The algorithm used in this code is referenced from the following: S
%   Gillijns et. al., "What Is the Ensemble Kalman Filter and How Well Does
%   it Work?" Proceedings of the 2006 American Control Conference,
%   Minneapolis, Minnesota, USA, June 14-16, 2006, pp 4448-4453.
%
%   ===================================================================   %


function [x_EnKF, x_tr] = OF_parallel_ensemblekfilter(T0,N,C,w,v,q,y_meas,sigma,tPlot,caseFolder_OF,solverName,tEnd,dt,scriptsFolder,solverRuns,solverCall) 
% OUTPUTS
%   x_EnKF: Predicted state forward in time (k+1) (with noise)
%
%   x_tr: Truth state (with noise)
%
%   MSE: mean-squred error between x_tr and x_EnKF for error quantification
%
% INPUTS
%   T0:     Initial temperature at t=0
%
%    N:     Number of cells
%
%    C:     Observation operator matrix maps the prior state vector into the 
%           vector of measurements (observables)---> y_k = C_k * x_k (nxn
%           matrix). If mapping the
%           state, then C = eye(nxn).
%    w:     Standard deviation of model noise with a Gaussian distribution
%
%    v:     Standard deviation of measurement noise with a Gaussian
%           distribution
%
%    q:     Number of ensemble members (greater q, less noise in x_EnKF
%           output but greater computational requirements)
%
%  y_meas:  Measurement matrix (n,k) ---> rows are number of states (or
%           cells), columns are number of time iterations (k). If no
%           measurement is available, use "nan"). Initiate as y_meas = nan
%           * ones(N,tEnd) then fill in availble measurements at
%           corresponding times.
%
%   sigma:  standard deviation of random sampling error noise to initialize
%           the ensemble
%
%   tplot:  time interval "k" to be plotted to show x_EnKF, x_tr and meas
%
%   caseFolder_OF: directory path to the OpenFOAM case folder
%                   i.e.) caseFolder_OF = '~/OpenFOAM/1D_Heat_Conduction/'
%
%   solverName: name of the solver to be called in OpenFOAM
%               i.e.) solverName = 'myLaplacianFoam'
%
%   tEnd:   end time
%
%   dt:     time step 
%
%   scriptsFolder: directory path to the MATLAB scripts
%                   i.e.) scriptsFolder = '~/dataAssimilation'
%
%   solverRuns: number of times OF solver should run per EnKF time
%               advancement
%   ===================================================================   %

% ======================== INITIALIZATION =============================   %
num_iterations = tEnd/(dt*solverRuns);
[y_rows,y_col] = size(y_meas);     % # rows, # col (nxk)
[C_rows,C_col] = size(C);          % # rows (m), # col (nx1)
x_est = T0 .* ones(N,q);       % generate ensemble for estimate state
if y_col ~= num_iterations         % y_meas must have a measurement or "nan" for each state (cell) and time iteration
    error('Must have a measurement for each individual time iteration. y_meas must have a column for each time iteration, k. Use "nan" if no measurement available.')
end

%   ==================== 1.ITERATE THROUGH TIME =======================   %

for k = 2:num_iterations + 1             % time loop (k=1 --> t=0, k=2 --> t=1*dt, k=3 --> t=2*dt)
   
    
   %  =============== 2.ITERATE THROUGH ENSEMBLE MEMBERS ==============   %
   for j = 1:q                              % loop through ensemble members one at a time
     W(:,j) = w .* randn(N,1);         % random Gaussian distributed noise with standard deviation input "w" (array size nx1)
     V(:,j) = v .* randn(C_rows,1);         % random Gaussian distributed noise with standard deviation input "z" (array size nx1)
     samp_err(:,j) = sigma .* randn(N,1);     % sample error used to generate initial ensemble below
     
     
     % ======== INITIALIZE THE ENSEMBLE WITH RANDOM SAMPLE ERROR ======== %
     if k == 2          % sample error only added during first time iteraiton
        x_est(:,j) = x_est(:,j) + samp_err(:,j);    % add sample error to each ensemble member
     end
     
     y_for(:,j) = C * x_est(:,j);           % forecast measurement (nxn * nxq = nxq matrix) 
     y(:,j) = y_meas(:,k-1) + V(:,j);         % add noise to measurements y_k,i = y_k + v_k,i (equations 3.6) (nxq matrix)
   end
   
   % Remove noise from B.C. cells (1st and last row in A matrix are B.C.s)
   W(1,:) = 0;
   W(end,:) = 0;
   V(1,:) = 0;
   V(end,:) = 0;
     
   % Replace "NaN" in y_meas matrix with the forecasted measurement, y_for
   % (when no measurement is available)
   if sum( isnan(y_meas(:,k-1)) ) >= 1              % if y_meas has one or more values that are NaN, then
   [row,col] = find(isnan(y_meas(:,k-1)));          % find row,col index of NaN elements
     for index = 1:length(row)
        y(row(index),:) = y_for(row(index),:);   % replace values of NaN from y_meas with y_for forecasted values so that these measurements are not assimilated with the K matrix
     end
   end
   
   
   %  =============== 3.AVERAGE ALL ENSEMBLE MEMBERS ==================   %
   x_estbar = mean(x_est,2);                % mean of ensemble of forecasted state (nx1)  
   y_forbar = mean(y_for,2);                % mean of ensemble of forecasted measurement (nx1)  
   
   
   %  ================= 4.ENSEMBLE ERROR MATRICES =====================   %
   for j = 1:N
     Ex(j,:) = [x_est(j,:) - x_estbar(j)];  % ensemble error matrix --> (n x q) matrix
   end 
   
   for j = 1:C_rows
     Ey(j,:) = [y_for(j,:) - y_forbar(j)];  % ensemble of output error matrix --> (p x q) matrix
   end
   
   %  ================= 5.ESTIMATE COVARIANCE MATRICES ================   %
   Pxy = Ex*Ey'/(q-1);                      % covariance matrix (nxq * qxp = nxp, q is # ensembles, p is # measurements, n is # of states)
   Pyy = Ey*Ey'/(q-1);                      % covariance matrix (pxq * qxp = pxp, q is # ensembles, p is # measurements, n is # of states)
   
   %  ==================== 6.KALMAN GAIN MATRIX =======================   %
   % (Pyy)^-1 is found using a pseudo-SVD decomposition to prevent
   % singularity numerical error
%    PyyInv = pinv(Pyy);
%    K = Pxy * PyyInv;
    K = eye(N);
   
   %  ======================= 7.ANALYSIS STEP =========================   %
   x_est = x_est + K * (y - y_for);         % new state estimate (nxp * pxq = nxq)
   
   % Send analsis step output to OpenFOAM and run the solver for one dt
   cd(scriptsFolder);   % change directory to folder with Matlab scripts
   tFolder = write_controlDict(k,caseFolder_OF,solverName,dt,tEnd,solverRuns); % write new controlDict file to start at the current time t
   varname = 'T';
   T = zeros(N,q);    % rows are each second of time, columns are cell centers
   for j = 1:q   % iterate through for each ensemble member
       changeFolderOF(varname,tFolder,caseFolder_OF);      % change directories to current OF time folder to create 'T' file from analysis step
       Matlab2OF(x_est,j,varname,tFolder,scriptsFolder);   % generate the T file for OpenFOAM from analysis step
       Tsolver = Matlab_callOF_myLaplacianFoam(k,varname,caseFolder_OF,solverCall,dt,solverRuns);  % change to t+dt time folder, run the solver for one time step
       T(:,j) = Tsolver;    % store values from OF solver output to a matrix (each row is a new cell, each column a new ensemble)
   
   %  ====================== 8. FORECAST STEP =========================   %
                              
     x_est(:,j) = T(:,j) + W(:,j);  % forecast state forward in time (nxn * nxq = nxq)
   end
   
   %  ====================== OUTPUT VARIABLES =========================   %
    x_EnKF(:,k) = mean(x_est,2);             % EnKF state output (each row is a state (cell), each column is an iteraiton in time)
    x_tr(:,k) = mean(T,2);     % true state = OpenFOAM output without noise
    x_EnKF(:,1) = T0;                 % set the initial condition for k=1
    x_tr(:,1) = T0;
end
% ======================== POST-PROCESSING ============================   %
ensemblekfilter_plots(tPlot,q,w,v,y_meas,x_EnKF,x_tr,dt,solverRuns)
end
