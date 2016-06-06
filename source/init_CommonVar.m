function [dirInfo,genericInfo]=init_CommonVar(fileName)%, callTypeTag)

fid=fopen(fileName);
C1=textscan(fid,'%s');
% variables
dirInfo.pathEMI=cell2mat(C1{1}(1));          % path for emission files
dirInfo.pathFOD=cell2mat(C1{1}(2));          % path for coordinate, policy application domain, AQI computation domain, population
dirInfo.pathFOO=cell2mat(C1{1}(3));          % path for optimization configuration file
dirInfo.pathOCM=cell2mat(C1{1}(4));          % path for measures DB
dirInfo.pathLST=cell2mat(C1{1}(5));          % path for file with techs list (to create NH3 constraints)
dirInfo.pathTRD=cell2mat(C1{1}(6));          % path for traffic duplication
dirInfo.pathDEF=cell2mat(C1{1}(7));          % path for AQI definition
dirInfo.pathPM =cell2mat(C1{1}(8));          % PM10: relation to convert "yearly average" to "# daily threshold exceedances"
dirInfo.pathOUF=cell2mat(C1{1}(9));          % path for output file
dirInfo.pathAR=cell2mat(C1{1}(10));          % path for areal_ratio variable (used for ASPA, for cells close to domain boundary)
dirInfo.pathADS=cell2mat(C1{1}(11));         % path for scenario analysis file, containing ms-pollutant % emission reductions

%[pathstr,name,ext] = fileparts(dirInfo.pathFOO); 
%dirInfo.pathCONF=strcat(pathstr, filesep, 'aggregated_configuration_', callTypeTag, '.txt');         % path for scenario analysis file, containing ms-pollutant % emission reductions

nocID=cell2mat(C1{1}(12));           % no control technology ID
genericInfo.AQINum=str2num(cell2mat(C1{1}(13))); % number of AQIs to consider
genericInfo.nocID=str2num(nocID);                % nocID to be transformed in number

fclose(fid);

end

