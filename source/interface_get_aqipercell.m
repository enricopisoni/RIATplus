function [ aqipercell ] = interface_get_aqipercell(input_rete2, NN, refInfo, commonDataInfo, periodIndex, aqiIndex)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_aqipercellFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'input_rete2,');
 stringToRun=strcat(stringToRun,'NN,');
 stringToRun=strcat(stringToRun,'refInfo,');
 stringToRun=strcat(stringToRun,'commonDataInfo,');
 stringToRun=strcat(stringToRun,'periodIndex,');
 stringToRun=strcat(stringToRun,'aqiIndex');
 stringToRun=strcat(stringToRun,')');
 [aqipercell]=eval(stringToRun);

end