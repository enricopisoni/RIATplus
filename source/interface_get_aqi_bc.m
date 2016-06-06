function [ bc ] = interface_get_aqi_bc(mathInfo,index)

stringToRun='';
stringToRun=strcat(stringToRun, mathInfo.get_aqi_bcFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'mathInfo,');
stringToRun=strcat(stringToRun,'index');
stringToRun=strcat(stringToRun,')');
[bc]=eval(stringToRun);

end