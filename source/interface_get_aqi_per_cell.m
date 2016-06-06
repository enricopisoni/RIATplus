function [ aqipercell ] = interface_get_aqi_per_cell(mathInfo,index1,index2)

stringToRun='';
stringToRun=strcat(stringToRun, mathInfo.get_aqi_per_cellFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'mathInfo,');
stringToRun=strcat(stringToRun,'index1,');
stringToRun=strcat(stringToRun,'index2');
stringToRun=strcat(stringToRun,')');
[aqipercell]=eval(stringToRun);

end