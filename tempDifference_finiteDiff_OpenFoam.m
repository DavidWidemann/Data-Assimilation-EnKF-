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

%% =================  FINITE DIFFERENCE METHOD  ======================= %%
% Here, a transient, 1D metal rod is discretized using the finite element 
% method into N cells. The transient heat conduction
% equation was spatially discretized using a central diffference method and
% temporally discretized using a first-order explicit Euler scheme. An
% absolute valued sinusoidal source term with a time dependence is utilized
% for heat addition at locations L/3 and L/6 of the rod. Note, T(0,t) =
% T(L,t) = T(x,0) = 300 K (BCs = IC = 300K). This temperature distrubtion
% along the bar as a function of time is used as experimental measurements
% for assimilation into the EnKF below.


%---------Boundary & Initial Conditions, Matrix Initialization---------%
T = zeros(length(t),N+1);   % each row is new time step, each column is new cell edge
Ti = 300;
T(:,1) = Ti;        % boundary condition
T(:,end) = Ti;      % boundary condition
T(1,:) = Ti;        % Initial condition

B = zeros(N+1,1);    % souce input mapping to cells 0.33L and 0.67L
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
%% ===========================  OPENFOAM  ============================= %%
% Here, Matlab calls OpenFOAM and runs a modified Laplacian solver called
% (myLaplacian). This file specifically just imports data from the folder
% containing the case and reads the temperatures from OpenFOAM for t_ini to
% t_end to a matrix for comparision against temperatures obtained from the
% finite difference method implemented above. It was found that the finite
% volume method used in OpenFOAM to solve the transient heat conduction
% equation yields temperatures most nearly the same as the a spatial
% central difference scheme with the maximum tempearture difference being
% 0.8 K, respectively.


scriptsFolder = '~/dataAssimilation';
caseFolder = '~/OpenFOAM/jns14008-5.x/1D_heatConductionEnKF';
cd(scriptsFolder);
fontsize = 17;
nheader = 22;
cd(caseFolder);

T2 = zeros(tEnd+1,N);    % rows are each second of time, columns are cell centers
for time = 1:tEnd
    str2 = sprintf('%d',time);
    cd(str2);
    fid = fopen('T');
    temp = textscan(fid,'%f', 'headerlines',nheader);
    fclose(fid);
    T2(time+1,:) = temp{1};
    cd ../;
end
T2(1,:) = 300;       % initialize temperture [K] at t=0, 0 <= x <= L 
[ghost,numCells] = size(T2);
cellCenterIDs = 1:length(numCells);

% PLOT OF TEMPERATURE DISTRIBUTION VS TIME
figure();
imagesc(cellCenterIDs,t,T2);
xlabel('Cell Center ID')
ylabel('Time (s)')
c = colorbar;
c.Label.String = 'Temperature (K)';
str3 = sprintf('1D Transient Heat Conduction in a Rod: q_{source}= 5e+07 * sin(t) W/m^3');
title(str3);

%% ======================  PLOT T DIFFERERNCE  ======================== %%
% 
% T1 = T(:,1:end-1);
% T1 = T1(1:1/dt:length(T1),:); 
% Tdiff = abs( T1 - T2 );
% 
% figure();
% imagesc(cellCenterIDs,t,Tdiff)
% xlabel('Cell Center ID')
% ylabel('Time (s)')
% c = colorbar;
% c.Label.String = 'Temperature (K)';
% title('Absolute Temperature Difference');
%% Test OpenFOAM Finite Volume Order of Convergence
T1 = T(1:1/dt:length(T),:); 
T1 = (T1(:,1:end-1) + T1(:,2:end) )/2;


% Use Taylor series expansion to calcualte dT/dx
for k = 1:length(T1)
    for i = 1:N
        T1(k,i) = ( T(k,i+1)-T(k,i) )/dx; % Matlab finite difference temp.
        T2(k,i) = ( T(k,i+1)-T(k,i) )/dx; % OpenFOAM finite volume temp.
    end
end

% Use Taylor series expansion again to calcualte d^2T/dx^2
for k = 1:length(T1)
    for i = 1:N
        T1(k,i) = ( T(k,i+1)-T(k,i) )/dx; % Matlab finite difference temp.
        T2(k,i) = ( T(k,i+1)-T(k,i) )/dx; % OpenFOAM finite volume temp.
    end
end

figure();
hold on;
plot(1:N,T1(end,:),':');
plot(1:N,T2(end,:),'--');
legend('finite difference','OpenFOAM finite volume'); legend box off;
hold off;

% This plot shows that the central finite difference and 2nd order finite volume
% schemes are exactly the same for 1D cases















