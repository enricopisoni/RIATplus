function [commonDataInfo]=interface_new_fill_common_data_info(commonDataInfo, periodIndex, aqiIndex, aggregationInfo)

if isequal(aggregationInfo.type, 'FIRSTGUESS')
    commonDataInfo=firstguess_fill_commonDataInfo(commonDataInfo, periodIndex, aqiIndex);
else
    commonDataInfo=quadrant_fill_commonDataInfo(commonDataInfo, periodIndex, aqiIndex);
end

end
