function [ optSet ] = interface_get_aqi_nnOptSet(refInfo,index)

stringToRun='';
stringToRun=strcat(stringToRun, refInfo.get_aqi_nnOptSetFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'refInfo,');
stringToRun=strcat(stringToRun,'index');
stringToRun=strcat(stringToRun,')');
[optSet]=eval(stringToRun);

end