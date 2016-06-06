function [ D ] = interface_get_DSuperSet(refInfo)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_DSuperSetFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'refInfo');
 stringToRun=strcat(stringToRun,')');
 [D]=eval(stringToRun);

end