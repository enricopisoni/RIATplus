function [ resEmis ] = interface_rebuildOrderedEmis(emis, refInfo, commonDataInfo)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_rebuildOrderedEmisFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'emis,');
 stringToRun=strcat(stringToRun,'refInfo,');
 stringToRun=strcat(stringToRun,'commonDataInfo');
 stringToRun=strcat(stringToRun,')');
 [resEmis]=eval(stringToRun);

end