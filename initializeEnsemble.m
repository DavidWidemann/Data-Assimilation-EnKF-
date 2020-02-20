function [x_ini_sootFv] = initializeEnsemble(varname,startTime,caseFolder_OF)

    %% ------------- Varname1 ------------- %%
    % # of header lines to skip
    nheader = 22;           % # of header lines to skip in OF file(s)
    nheaderBC1 = 22313;     % 1st boundary condition values
    nheaderBC2 = 22334;
    nheaderBC3 = 22423;
    nheaderBC4 = 22497;
    nheaderBC5 = 22667;
    nheader = [nheader nheaderBC1 nheaderBC2 nheaderBC3 nheaderBC4 nheaderBC5];

    % Change to start time folder
    str = sprintf('%s/%.15g',caseFolder_OF,startTime);
    cd(str);
    
    % Read OpenFOAM output into a Matlab array/matrix
    fid = fopen(varname1);
        data1 = textscan(fid,'%f', 'headerlines',nheader(1));
    fclose(fid);
    
    fid = fopen(varname1);
        data2 = textscan(fid,'%f', 'headerlines',nheader(2));
    fclose(fid);
    
    fid = fopen(varname1);
        data3 = textscan(fid,'%f', 'headerlines',nheader(3));
    fclose(fid);
    
    fid = fopen(varname1);
        data4 = textscan(fid,'%f', 'headerlines',nheader(4));
    fclose(fid);
    
    fid = fopen(varname1);
        data5 = textscan(fid,'%f', 'headerlines',nheader(5));
    fclose(fid);
    
    fid = fopen(varname1);
        data6 = textscan(fid,'%f', 'headerlines',nheader(6));
    fclose(fid);
    
    x_ini_sootFv = [data1{1}; data2{1}; data3{1}; data4{1}; data5{1}; data6{1}] ;
    
end