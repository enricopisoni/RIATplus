function [ aqipercell ] = do_get_aqi_per_cell(input_rete2, refInfo, periodIndex, aqiIndex)

 stringToRun='';
 stringToRun=strcat(stringToRun, refInfo.get_do_aqi_per_cellFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'emissioni,');
 stringToRun=strcat(stringToRun,'refInfo,');
 stringToRun=strcat(stringToRun,'periodIndex,');
 stringToRun=strcat(stringToRun,'aqiIndex');
 stringToRun=strcat(stringToRun,')');
 [aqipercell]=eval(stringToRun);

end