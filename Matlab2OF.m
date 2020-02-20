
%% Script to take data from MATLAB and convert to OpenFOAM
%% Xinyu Zhao 2017, Joseph Squeo 2019
function Matlab2OF(x_est,j,varname,tFolder,scriptsFolder,x0)
%% Write data in OF format
NDATA = size(x_est,1);
Tdata = x_est(:,j);

    %fprintf('%s - q%.15g\n',varname,j);
    fid = fopen(varname,'w'); % write permission
    fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\\n');
    fprintf(fid,'| =========                 |                                                 |\n');
    fprintf(fid,'| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n');
    fprintf(fid,'|  \\\\    /   O peration     | Version:  5.x                                   |\n');
    fprintf(fid,'|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n');
    fprintf(fid,'|    \\\\/     M anipulation  |                                                 |\n');
    fprintf(fid,'\\*---------------------------------------------------------------------------*/\n');
    fprintf(fid,'FoamFile\n');
    fprintf(fid,'{\n');
    fprintf(fid,'    version     2.0;\n');
    fprintf(fid,'    format      ascii;\n');
    if strcmp(varname,'U')
        fprintf(fid,'    class       volVectorField;\n');
    else
        fprintf(fid,'    class       volScalarField;\n');
    end
    fprintf(fid,'    location    "%.15g";\n',tFolder);
    fprintf(fid,'    object       %s;\n',varname);
    fprintf(fid,'}\n');
    fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');
    fprintf(fid,'\n');
    %     fprintf(fid,'internalField   uniform %.15g;\n',Yox(k));
    if strcmp(varname,'p')
        fprintf(fid,'dimensions      [1 -1 -2 0 0 0 0];\n');
        fprintf(fid,'\n');
        fprintf(fid,'internalField   nonuniform List<scalar>\n');
        fprintf(fid,'%d\n',NDATA);
        fprintf(fid,'(\n');
        fprintf(fid,'%.2f\n', N(:,iv-4));
        fprintf(fid,');\n');
    elseif strcmp(varname,'T')
        fprintf(fid,'dimensions      [0 0 0 1 0 0 0];\n');
        fprintf(fid,'\n');
        fprintf(fid,'internalField   nonuniform List<scalar>\n');
        fprintf(fid,'%d\n',NDATA);
        fprintf(fid,'(\n');
        fprintf(fid,'%.2f\n',Tdata);
        fprintf(fid,');\n');
    elseif strcmp(varname,'N2')
        fprintf(fid,'dimensions      [0 0 0 0 0 0 0];\n');
        fprintf(fid,'\n');
        fprintf(fid,'internalField   nonuniform List<scalar>\n');
        fprintf(fid,'%d\n',NDATA);
        fprintf(fid,'(\n');
        fprintf(fid,'%.5e\n',N(:,length(varnames)-4));
        fprintf(fid,');\n');
    else
        fprintf(fid,'dimensions      [0 0 0 0 0 0 0];\n');
        fprintf(fid,'\n');
        fprintf(fid,'internalField   nonuniform List<scalar>\n');
        fprintf(fid,'%d\n',NDATA);
        fprintf(fid,'(\n');
        fprintf(fid,'%.5e\n',N(:,iv-4));
        fprintf(fid,');\n');
    end 
    fprintf(fid,'\nboundaryField\n');
    fprintf(fid,'{\n');
    fprintf(fid,'    inlet\n');
    fprintf(fid,'    {\n');
       if strcmp(varname,'T')
        fprintf(fid,'        type            fixedValue;\n');
        fprintf(fid,'        value           uniform %.15g;\n',x0);
       else
            fprintf(fid,'    type           zeroGradient;\n');
       end
    fprintf(fid,'    }\n');
    fprintf(fid,'    exit\n');
    fprintf(fid,'    {\n');
    if strcmp(varname,'T')
        fprintf(fid,'        type            fixedValue;\n');
        fprintf(fid,'        value           uniform %.15g;\n',x0);
       else
            fprintf(fid,'    type           zeroGradient;\n');
    end
    fprintf(fid,'    }\n');
     fprintf(fid,'    frontAndBack\n');
    fprintf(fid,'    {\n');
    if strcmp(varname,'T')
        fprintf(fid,'        type            empty;\n');
        %fprintf(fid,'        value           uniform %.15g;\n',T0);
       else
            fprintf(fid,'    type           zeroGradient;\n');
    end
    fprintf(fid,'    }\n');
     fprintf(fid,'    topAndBottom\n');
    fprintf(fid,'    {\n');
    if strcmp(varname,'T')
        fprintf(fid,'        type            empty;\n');
        %fprintf(fid,'        value           uniform %.15g;\n',T0);
       else
            fprintf(fid,'    type           zeroGradient;\n');
    end
    fprintf(fid,'    }\n');
    fprintf(fid,'}\n');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'// ************************************************************************* //\n');
    
    fclose(fid);  
    cd(scriptsFolder);
end
