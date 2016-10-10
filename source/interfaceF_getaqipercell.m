function [aqi_per_cell]=interfaceF_getaqipercell(type, emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex)

if isequal(type, 'FIRSTGUESS')
    aqi_per_cell=firstguess_do_aqi_per_cell(emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex);
else
    aqi_per_cell=do_interface_get_aqi_per_cell(emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex);
end

end