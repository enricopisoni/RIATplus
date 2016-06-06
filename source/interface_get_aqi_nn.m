function [ nn ] = interface_get_aqi_nn(mathInfo,index)

stringToRun='';
stringToRun=strcat(stringToRun, mathInfo.get_aqi_nnFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'mathInfo,');
stringToRun=strcat(stringToRun,'index');
stringToRun=strcat(stringToRun,')');
[bc]=eval(stringToRun);

end