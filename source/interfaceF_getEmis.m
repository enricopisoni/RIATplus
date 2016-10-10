function [emissioni]=interfaceF_getEmis(type, emis, emisP, NN, bcSet, aggregationInfo, emisType, periodIndex, aqiIndex)

if isequal(type, 'FIRSTGUESS')
    bcSet=0;
    NN=0;
    aggregationInfo.mathInfo=0;
    % 20160420 fg version
    emissioni=firstguess_buildEmission(emis, emisP, NN, bcSet, aggregationInfo.mathInfo, aggregationInfo,emisType,periodIndex,aqiIndex);
else
    emissioni=interface_buildEmission(emis, emisP, bcSet, NN, aggregationInfo.mathInfo,emisType,periodIndex,aqiIndex);
end

end