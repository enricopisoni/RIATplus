function [ nn ] = interface_get_nnSuperSet_indexed(refInfo, periodIndex, aqiIndex)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_nnSuperSetIndexedFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'refInfo');
 stringToRun=strcat(stringToRun,')');
 [nn]=eval(stringToRun);

end