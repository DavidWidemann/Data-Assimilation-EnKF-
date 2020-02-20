
%% Calling OpenFOAM from Matlab:
function Tsolver = Matlab_callOF_myLaplacianFoam(k,varname,caseFolder_OF,solverCall,dt,solverRuns)
% INPUTS:
%           T: temperature matrix
%           q: number of ensembles 
%        tEnd: end time
%           N: number of cells (or states)
    
    % INITIALIZATION
    nheader = 22;       % # of header lines to skip in OF file(s)
    %callSolver = sprintf('%s',solverName);    % Solver name in OpenFOAM
    cd(caseFolder_OF);      % change directory to the case foler
    
    % CALL SOLVER IN OPENFOAM
    [x,y] = unix(solverCall);   % call and run the solver
    
    % CHANGE TO THE NEW TIME FOLDER TO STORE SOLVER OUTPUT
    tFolderNew = solverRuns*dt*(k-1);   % (k-1) because want new time folder
    if isreal(tFolderNew) && tFolderNew >= 0
        str1 = sprintf('%s/%.15g',caseFolder_OF,tFolderNew);
        cd(str1);
    else
        error('Cannot change to specified time folder!\n Time Folder is not a real number, a negative value or not recognized.\n');
    end
    
%     if parallel == 1
%         reconstruct = sprintf('reconstructPar -times ''%g:%g''',tFolder,tFolderNew);
%         unix(reconstruct);
%     end
    
    % READ SOLVER OUTPUT(S) INTO MATLAB AS MATRICES/VECTORS
    fid = fopen(varname);
    temp = textscan(fid,'%f', 'headerlines',nheader);
    fclose(fid);
    Tsolver(:,1) = temp{1};
 
end
