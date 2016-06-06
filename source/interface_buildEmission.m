function [ emis ] = interface_buildEmission(emis,emisP,NN,bcSet,refInfo,aggInfo,refData,emissionType,periodIndex,aqiIndex,isAggregated,isDelta)
%function [ emis ] = interface_buildEmission(emis,emisP,NN,bcSet,refInfo,aggInfo,emissionType,periodIndex,aqiIndex,isAggregated,isDelta)

 stringToRun='';
 stringToRun=strcat(stringToRun, refData.get_buildEmissionFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'emis,');
 stringToRun=strcat(stringToRun,'emisP,');
 stringToRun=strcat(stringToRun,'NN,');
 stringToRun=strcat(stringToRun,'bcSet,');
 stringToRun=strcat(stringToRun,'refInfo,');
 stringToRun=strcat(stringToRun,'aggInfo,');
 stringToRun=strcat(stringToRun,'emissionType,');
 stringToRun=strcat(stringToRun,'periodIndex,');
 stringToRun=strcat(stringToRun,'aqiIndex');
 stringToRun=strcat(stringToRun,')');
 [emis]=eval(stringToRun);

end