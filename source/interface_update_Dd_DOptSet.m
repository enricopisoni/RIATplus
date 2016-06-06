function [res] = interface_update_Dd_DOptSet(D, d, refInfo, index1)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.update_Dd_DOptSetFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'D, ');
 stringToRun=strcat(stringToRun,'d, ');
 stringToRun=strcat(stringToRun,'refInfo,');
 stringToRun=strcat(stringToRun,'index1');
 stringToRun=strcat(stringToRun,')');
 res=eval(stringToRun);

end