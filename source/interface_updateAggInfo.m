function [ ncommonDataInfo, naggregationInfo  ] = interface_fillcomputesol(commonDataInfo,refInfo,aggregationInfo,periodIndex,aqiIndex)

stringToRun='';
stringToRun=strcat(stringToRun, refInfo.get_fillcomputesolFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'commonDataInfo,');
stringToRun=strcat(stringToRun,'refInfo,');
stringToRun=strcat(stringToRun,'aggregationInfo,');
stringToRun=strcat(stringToRun,'periodIndex,');
stringToRun=strcat(stringToRun,'aqiIndex');
stringToRun=strcat(stringToRun,')');
[ncommonDataInfo, naggregationInfo]=eval(stringToRun);

end