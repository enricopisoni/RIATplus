function [flag] = interface_is6Input(refInfo,periodIndex,aqiIndex)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.is6InputFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'refInfo,');
 stringToRun=strcat(stringToRun,'periodIndex,');
 stringToRun=strcat(stringToRun,'aqiIndex');
 stringToRun=strcat(stringToRun,')');
 flag=eval(stringToRun);

end