function [emis]=interfaceF_rebuildOrderedEmis(type, refEp, geometryIntermediateData, commonDataInfo)

if isequal(type, 'FIRSTGUESS')
    emis=FG_rebuildOrderedEmis(refEp, geometryIntermediateData, commonDataInfo);
else
    % 20160420 nn/quadrant version
    emis=interface_rebuildOrderedEmis(refEp, geometryIntermediateData, commonDataInfo);
end

end