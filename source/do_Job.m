function [intermediateResult, resultInfo]=do_Job(refDataInfo, commonData, precursor, indicators, x, y, nx, ny, ...
    precName, mathIndex, splitResult)
% function [intermediateResult, resultInfo]=do_Job(refDataInfo, commonData, precursor, indicators, x, y, nx, ny, totalCells, ...
%     optimizerValues, optimizerCondition)

 %intermediateResult=interface_prepare(refDataInfo.prepareInfo, commonData, precursor, indicators, x, y, nx, ny, totalCells, optimizerValues, optimizerCondition); 
 intermediateResult=interface_prepare(refDataInfo, commonData); 
 resultInfo=interface_finalize(refDataInfo, intermediateResult, commonData, precursor, indicators, x, y, nx, ny, precName, mathIndex, splitResult);
 %resultInfo=interface_finalize(refDataInfo, intermediateResult, commonData, precursor, indicators, x, y, nx, ny, totalCells, optimizerValues, optimizerCondition);

end