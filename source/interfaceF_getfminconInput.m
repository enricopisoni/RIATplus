function [D, bc, nn]=interfaceF_getfminconInput(type, aggregationInfo, ii)

if isequal(aggregationInfo.type, 'FIRSTGUESS')
    D=firstguess_get_aqi_D(aggregationInfo.geometryIntermediateData, ii);
    %                 sum(sum(D.D))
    %                 sum(D.d)
    bc=1;
    % 20160421 First guess setting
    nn=1;
else
    % 20160421 quadrant/NN setting
    D=interface_get_aqi_D(aggregationInfo.geometryIntermediateData, ii);
    %cyclebc, aggInfo, commonDataInfo, ii, jj
    %bc: delta or abs values
    % 20160421 First guess setting
    %bc: delta or abs values
    % 20160421 quadrant/NN setting
    bc=interface_get_aqi_bc(aggregationInfo.mathIntermediateData, ii);
    % 20160421 quadrant/NN setting
    nn=interface_get_aqi_nnOptSet(aggregationInfo.mathIntermediateData, ii);
end

end