function [flag] = interface_isAggregatedEmis(refInfo, periodIndex, aqiIndex)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.isAggregatedEmisFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'refInfo,');
 stringToRun=strcat(stringToRun,'periodIndex,');
 stringToRun=strcat(stringToRun,'aqiIndex');
 stringToRun=strcat(stringToRun,')');
 flag=eval(stringToRun);

end