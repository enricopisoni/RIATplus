function [ D,d ] = interface_get_Dd_DOptSet(refInfo,index1)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_Dd_DOptSetFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'refInfo,');
 stringToRun=strcat(stringToRun,'index1');
 stringToRun=strcat(stringToRun,')');
 [D,d]=eval(stringToRun);

end