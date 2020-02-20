function BC_Data = XC2H2BoundaryData(caseFolder_OF,startTime,varname)

str1 = sprintf('%s/%.15g',caseFolder_OF,startTime);
cd(str1);

% READ SOLVER OUTPUT(S) INTO MATLAB AS MATRICES/VECTORS
    fid = fopen(varname);
        data1 = textscan(fid,'%f', 'headerlines',22313);
    fclose(fid);
    
    fid = fopen(varname);
        data2 = textscan(fid,'%f', 'headerlines',22334);
    fclose(fid);
    
    fid = fopen(varname);
        data3 = textscan(fid,'%f', 'headerlines',22428);
    fclose(fid);
    
    fid = fopen(varname);
        data4 = textscan(fid,'%f', 'headerlines',22598);
    fclose(fid);
    
    BC_Data = [data1{1}; data2{1}; data3{1}; data4{1}] ;