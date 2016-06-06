function [intermediateResult, resultInfo]=firstguess_do_Job(refDataInfo, commonData, precursor, indicators, x, y, nx, ny, ...
    precName, mathIndex, splitResult)
% function [intermediateResult, resultInfo]=do_Job(refDataInfo, commonData, precursor, indicators, x, y, nx, ny, totalCells, ...
%     optimizerValues, optimizerCondition)

 %intermediateResult=interface_prepare(refDataInfo.prepareInfo, commonData, precursor, indicators, x, y, nx, ny, totalCells, optimizerValues, optimizerCondition); 
 %intermediateResult=FG_interface_prepare(refDataInfo, commonData); 
 %resultInfo=FG_interface_finalize(refDataInfo, intermediateResult, commonData, precursor, indicators, x, y, nx, ny, precName, mathIndex, splitResult);
 %init...
 intermediateResult=firstguess_Prepare_init(refDataInfo, commonData); 
 resultInfo=firstguess_Finalize(refDataInfo, intermediateResult, commonData, precursor, indicators, x, y, nx, ny, precName, mathIndex, splitResult);
 %resultInfo=interface_finalize(refDataInfo, intermediateResult, commonData, precursor, indicators, x, y, nx, ny, totalCells, optimizerValues, optimizerCondition);

end