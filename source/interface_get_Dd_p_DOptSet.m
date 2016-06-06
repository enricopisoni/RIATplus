function [ Dp,dp ] = interface_get_Dd_p_DOptSet(refInfo,index1)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_Dd_p_DOptSetFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'refInfo,');
 stringToRun=strcat(stringToRun,'index1');
 stringToRun=strcat(stringToRun,')');
 [Dp,dp]=eval(stringToRun);

end