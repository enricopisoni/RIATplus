function [localnnSuperSet, localDSuperSet, localbcSuperSet]=interfaceF_getlocalSuperSet(type, aggregationInfo)

%if isequal(aggregationInfo.type, 'FIRSTGUESS')
if isequal(type, 'FIRSTGUESS')
    localnnSuperSet=0;
    localDSuperSet=firstguess_get_aqi_D(aggregationInfo.geometryIntermediateData, 1);
    localbcSuperSet=0;
else
    localnnSuperSet=interface_get_nnSuperSet(aggregationInfo.mathIntermediateData);
    localDSuperSet=interface_get_DSuperSet(aggregationInfo.geometryIntermediateData);
    localbcSuperSet=interface_get_bcSuperSet(aggregationInfo.mathIntermediateData);
end

end

