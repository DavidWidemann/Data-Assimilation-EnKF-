
%% Ccript to take data from MATLAB and convert to OpenFOAM
%% Xinyu Zhao 2017, Joseph Squeo 2019
function tFolder = write_controlDict(k,caseFolder_OF,solverName,dt,tEnd,solverRuns)
%% Change directories
tFolder = solverRuns*dt*(k-2);
str1 = sprintf('%s/system/',caseFolder_OF);
cd(str1);
maxCourant = 0.5;

%% Write controlDict in OF format

    fid = fopen('controlDict','w'); % write permission
    fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\\n');
    fprintf(fid,'| =========                 |                                                 |\n');
    fprintf(fid,'| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n');
    fprintf(fid,'|  \\\\    /   O peration     | Version:  2.3.0                                 |\n');
    fprintf(fid,'|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n');
    fprintf(fid,'|    \\\\/     M anipulation  |                                                 |\n');
    fprintf(fid,'\\*---------------------------------------------------------------------------*/\n');
    fprintf(fid,'FoamFile\n');
    fprintf(fid,'{\n');
    fprintf(fid,'    version     2.0;\n');
    fprintf(fid,'    format      ascii;\n');
    fprintf(fid,'    class       dictionary;\n');   
    fprintf(fid,'    location    "system";\n');
    fprintf(fid,'    object       controlDict;\n');
    fprintf(fid,'}\n');
    fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');
    fprintf(fid,'\n');

    fprintf(fid,'application     %s;\n',solverName);
    fprintf(fid,'\n');
    fprintf(fid,'startFrom       startTime;\n');
    fprintf(fid,'\n');
    fprintf(fid,'startTime       %.15g;\n',tFolder);
    fprintf(fid,'\n');
    if solverRuns == 1
        fprintf(fid,'stopAt          writeNow;\n');
    elseif solverRuns >= 2
        fprintf(fid,'stopAt          nextWrite;\n');

    else
        error('solverRuns must be a whole, positive number!');
    end
    fprintf(fid,'\n');
    fprintf(fid,'endTime         %.15g;\n',tEnd);
    fprintf(fid,'\n');
    fprintf(fid,'deltaT          %.15g;\n',dt);
    fprintf(fid,'\n');
    fprintf(fid,'writeControl    adjustableRunTime;\n');
    fprintf(fid,'\n');
    if solverRuns == 1
        fprintf(fid,'writeInterval   %.15g;\n',dt);
    elseif solverRuns >= 2
        fprintf(fid,'writeInterval   %.15g;\n',dt*solverRuns);
    else
        error('solverRuns must be a whole, positive number!');
    end
    fprintf(fid,'\n');
    fprintf(fid,'purgeWrite      0;\n');
    fprintf(fid,'\n');
    fprintf(fid,'writeFormat     ascii;\n');
    fprintf(fid,'\n');
    fprintf(fid,'writePrecision  10;\n');
    fprintf(fid,'\n');
    fprintf(fid,'writeCompression off;\n');
    fprintf(fid,'\n');
    fprintf(fid,'timeFormat      general;\n');
    fprintf(fid,'\n');
    fprintf(fid,'timePrecision   8;\n');
    fprintf(fid,'\n');
    fprintf(fid,'runTimeModifiable true;\n');
    fprintf(fid,'\n');
    fprintf(fid,'adjustTimeStep  no;\n');
    fprintf(fid,'\n');
    fprintf(fid,'maxCo           %.15g;\n',maxCourant);
    fprintf(fid,'\n');
    fprintf(fid,'functions\n');
    fprintf(fid,'{\n');
    fprintf(fid,'}\n');
    fprintf(fid,'\n');
    fprintf(fid,'// ************************************************************************* //\n');

    fclose(fid);  
    cd(caseFolder_OF);
end
