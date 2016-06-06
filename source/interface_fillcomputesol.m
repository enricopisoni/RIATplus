function [ ncommonDataInfo, naggregationInfo  ] = interface_fillcomputesol(commonDataInfo,aggregationInfo)

stringToRun='';
stringToRun=strcat(stringToRun, aggregationInfo.get_fillcomputesolFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'commonDataInfo,');
stringToRun=strcat(stringToRun,'aggregationInfo');
stringToRun=strcat(stringToRun,')');
[ncommonDataInfo, naggregationInfo]=eval(stringToRun);

end