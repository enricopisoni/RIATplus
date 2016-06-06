function [ bc ] = interface_get_bcSuperSet(refInfo)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_bcSuperSetFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'refInfo');
 stringToRun=strcat(stringToRun,')');
 [bc]=eval(stringToRun);

end