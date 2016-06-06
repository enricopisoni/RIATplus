%function [ aqi_per_cell ]= neuralNet_get_aqipercell(input_rete2, refData, commonDataInfo, periodIndex, aqiIndex)
function [ aqi_per_cell ]= firstguess_get_aqipercell(input, refData, commonDataInfo, periodIndex, aqiIndex)

 aqi_per_cell=input';
 %PATCH 20160601 EP
%  aqi_per_cell(isnan(aqi_per_cell))=[];

end
