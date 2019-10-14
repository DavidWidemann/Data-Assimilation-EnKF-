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
%   * VARIABLES *
%   k:      time iteration
%   A:      input model operator matrix
%   B:      input model source term mapping vector
%   uk:     input model source term values stored in an array
%   x_tr:   truth state
%   x_est:  analsysis and/or forecasted state
%   q:      # of ensemble members
%   C:      measurement mapping operator
%   x_ini:  inital condition of the state
%   w,v:    model and measurement Gaussian noise standard deviations
%   y_meas: measurement inputs of the state
%   num_iterations: number of time steps
%   sigma: standard deviaiton of random sample error for ensemble
%   K:      kalman gain matrix
%   y_for:  forecasted state measurements (y_for = C*x_est)
%   y:      measurements + Gaussian distributed random noise
%   x_estbar, y_estbar: average of ensembles
%   Ex, Ey: enseble error matrices for state and forecasted measurements
%   Pxx, Pyy: estimated covariance matrices
%   MSE:     mean-squared error
%   x_EnKF: predicted state forward in time from ensemble kalman filter

%   ===================================================================   %


function [x_EnKF, x_tr,MSE] = ensemblekfilter(A,C,x_ini,w,v,q,y_meas,num_iterations,B,uk,sigma,tPlot,dt) 

% OUTPUTS
%   x_EnKF: Predicted state forward in time (k+1) (with noise)
%
%   x_tr: Truth state (with noise)
%
%   MSE: mean-squred error between x_tr and x_EnKF for error quantification
%
% INPUTS
%   A:      Model operator matrix maps the previous state to the new state
%            ---> x_(k+1) = A_k * x_k (nxn matrix where n is size of
%            the state or number of cells in domain). This matrix would
%            contain the discretized governing equations.
%
%   C:      Observation operator matrix maps the prior state vector into the 
%           vector of measurements (observables)---> y_k = C_k * x_k (nxn
%           matrix). If mapping the
%           state, then C = eye(nxn).
%
%   x_ini:  Column vector of initial conditions of the state ---> input as
%           [1; 2; 1; 1...]
%           (each row represents a state (or cell). For a bar at T=300K and
%           100 cells, x_ini = 300*ones(100,1).
%
%   w:      Standard deviation of model noise with a Gaussian distribution
%
%   v:      Standard deviation of measurement noise with a Gaussian
%           distribution
%
%   q:      Number of ensemble members (greater q, less noise in x_EnKF
%           output but greater computational requirements)
%
%   y_meas: Measurement matrix (n,k) ---> rows are number of states (or
%           cells), columns are number of time iterations (k). If no
%           measurement is available, use "nan"). Initiate as y_meas = nan
%           * ones(N,tEnd) then fill in availble measurements at
%           corresponding times.
%
%   num_iterations: # of data assimilation cycles in time. If running for
%           100 seconds, then num_iterations = 100 or tEnd.
%
%   B:      input model source term mapping column vector to map source
%           term to specific locations ---> input as [0; 0; 1; 0 ...] for      
%           uk to affect the 3rd cell and so on
%
%   uk:     row vector containing source term values for each time
%           iteration sigma:  sample error standard deviation used to 
%           generate initial ensemble from x_ini
%
%   sigma:  standard deviation of random sample error to initialize
%           ensemble
%
%   tplot:  time t to be plotted 
%
%   dt:     time step 
%   ===================================================================   %


% ======================== INITIALIZATION =============================   %
[A_rows,A_col] = size(A);          % # rows, # col (nxn)
[y_rows,y_col] = size(y_meas);     % # rows, # col (nxk)
if A_rows ~= A_col
    error('A must be a square nxn matrix')
end
[C_rows,C_col] = size(C);          % # rows (m), # col (nx1)
x_tr = x_ini;                      % initialize true state at k=1
x_est = x_ini .* ones(A_rows,q);   % generate ensemble for estimate state

if y_col ~= num_iterations         % y_meas must have a measurement or "nan" for each state (cell) and time iteration
    error('Must have a measurement for each individual time iteration. y_meas must have a column for each time iteration, k. Use "nan" if no measurement available.')
end


%   ==================== 1.ITERATE THROUGH TIME =======================   %
for k = 2:num_iterations             % time loop (at k=1 , x_est = x_ini, the initial condition)
   
     x_tr(:,k) = A*x_tr(:,k-1) + B*uk(k-1) + w.*randn(A_rows,1); % compute truth state (nxn * nx1 = nx1)
   
    
   %  =============== 2.ITERATE THROUGH ENSEMBLE MEMBERS ==============   %
   for j = 1:q                              % loop through ensemble members one at a time
     W(:,j) = w .* randn(A_rows,1);         % random Gaussian distributed noise with standard deviation input "w" (array size nx1)
     V(:,j) = v .* randn(C_rows,1);         % random Gaussian distributed noise with standard deviation input "z" (array size nx1)
     samp_err(:,j) = sigma .* randn(A_rows,1);     % sample error used to generate initial ensemble below
     
     
     % ======== INITIALIZE THE ENSEMBLE WITH RANDOM SAMPLE ERROR ======== %
     if k == 2          % sample error only added during first time iteraiton
        x_est(:,j) = x_est(:,j) + samp_err(:,j);    % add sample error to each ensemble member
     end
     
     y_for(:,j) = C * x_est(:,j);           % forecast measurement (nxn * nxq = nxq matrix) 
     y(:,j) = y_meas(:,k) + V(:,j);         % add noise to measurements y_k,i = y_k + v_k,i (equations 3.6) (nxq matrix)
   end
   
   % Remove noise from B.C. cells (1st and last row in A matrix are B.C.s)
   W(1,:) = 0;
   W(end,:) = 0;
   V(1,:) = 0;
   V(end,:) = 0;
     
   % Replace "NaN" in y_meas matrix with the forecasted measurement, y_for
   % (when no measurement is available)
   if sum( isnan(y_meas(:,k)) ) >= 1              % if y_meas has one or more values that are NaN, then
   [row,col] = find(isnan(y_meas(:,k)));          % find row,col index of NaN elements
     for index = 1:length(row)
        y(row(index),:) = y_for(row(index),:);   % replace values of NaN from y_meas with y_for forecasted values so that these measurements are not assimilated with the K matrix
     end
   end
   
   
   %  =============== 3.AVERAGE ALL ENSEMBLE MEMBERS ==================   %
   x_estbar = mean(x_est,2);                % mean of ensemble of forecasted state (nx1)  
   y_forbar = mean(y_for,2);                % mean of ensemble of forecasted measurement (nx1)  
   
   
   %  ================= 4.ENSEMBLE ERROR MATRICES =====================   %
   for j = 1:A_rows
     Ex(j,:) = [x_est(j,:) - x_estbar(j)];  % ensemble error matrix --> (n x q) matrix
   end 
   
   for j = 1:C_rows
     Ey(j,:) = [y_for(j,:) - y_forbar(j)];  % ensemble of output error matrix --> (p x q) matrix
   end
   
   %  ================= 5.ESTIMATE COVARIANCE MATRICES ================   %
   Pxy = Ex*Ey'/(q-1);                      % covariance matrix (nxq * qxp = nxp, q is # ensembles, p is # measurements, n is # of states)
   Pyy = Ey*Ey'/(q-1);                      % covariance matrix (pxq * qxp = pxp, q is # ensembles, p is # measurements, n is # of states)

   
   %  ==================== 6.KALMAN GAIN MATRIX =======================   %
   % (Pyy)^-1 is found using a pseudo-SVD decomposition since Pyy is most nearly singular
   K = Pxy * pinv(Pyy);
   
   
   %  ======================= 7.ANALYSIS STEP =========================   %
   x_est = x_est + K * (y - y_for);         % new state estimate (nxp * pxq = nxq)
   
   
   %  ====================== 8. FORECAST STEP =========================   %
   for j = 1:q                              % number of ensemble members
     x_est(:,j) = A * x_est(:,j) + B * uk(k-1) + W(:,j);  % forecast state forward in time (nxn * nxq = nxq)
   end
   
   %  =============== s ==================   %
   MSE(k,:) = ( x_tr(:,k) - mean(x_est,2) ).^2 ; % mean squared error
   x_EnKF(:,k) = mean(x_est,2);             % EnKF state output (each row is a state (cell), each column is an iteraiton in time)
   x_EnKF(:,1) = x_ini;                 % set the initial condition for k=1
end


% ======================== POST-PROCESSING ============================   %

ensemblekfilter_plots(tPlot,q,w,v,y_meas,x_EnKF,x_tr,dt)
MSE_plots(MSE,q)
end
