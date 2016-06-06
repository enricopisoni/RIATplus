function [ nn ] = interface_get_nnOptSet(refInfo)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_nnOptSetFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'refInfo');
 stringToRun=strcat(stringToRun,')');
 [nn]=eval(stringToRun);

end