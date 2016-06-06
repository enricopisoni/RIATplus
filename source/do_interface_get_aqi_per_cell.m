function [aqi_per_cell] = do_interface_get_aqi_per_cell(emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex)

 stringToRun='';
 stringToRun=strcat(stringToRun, aggregationInfo.mathIntermediateData.get_aqi_per_cellFunction);
 stringToRun=strcat(stringToRun, '(');
 stringToRun=strcat(stringToRun,'emissioni,');
 stringToRun=strcat(stringToRun,'NN,');
 stringToRun=strcat(stringToRun,'aggregationInfo,');
 stringToRun=strcat(stringToRun,'commonDataInfo,');
 stringToRun=strcat(stringToRun,'periodIndex,');
 stringToRun=strcat(stringToRun,'aqiIndex');
 stringToRun=strcat(stringToRun,')');
 [aqi_per_cell]=eval(stringToRun);

end