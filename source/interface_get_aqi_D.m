function [ D ] = interface_get_aqi_D(mathInfo,index)

stringToRun='';
stringToRun=strcat(stringToRun, mathInfo.get_aqi_DFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'mathInfo,');
stringToRun=strcat(stringToRun,'index');
stringToRun=strcat(stringToRun,')');
[D]=eval(stringToRun);

end