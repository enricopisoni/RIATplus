function [mathIntermediateData]=interfaceF_prepareMathDataInfo(type, aggregationInfo, commonDataInfo)

if isequal(type, 'FIRSTGUESS')
    mathIntermediateData=0;
else
    mathIntermediateData=interface_prepare(aggregationInfo.mathDataInfo, commonDataInfo);
end

end

