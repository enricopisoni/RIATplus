function [ nn ] = interface_get_nnSuperSet(refInfo)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_nnSuperSetFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'refInfo');
 stringToRun=strcat(stringToRun,')');
 [nn]=eval(stringToRun);

end