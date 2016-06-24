function [intermediateResult resultInfo]=do_FullJob(refDataInfo, commonData, precursorValues, precursorNames, indicators, x, y, nx, ny, ...
    mathIndex, indoptrep, doSparse, doFilter)
% function [intermediateResult, resultInfo]=do_Job(refDataInfo, commonData, precursor, indicators, x, y, nx, ny, totalCells, ...
%     optimizerValues, optimizerCondition)

 %intermediateResult=interface_prepare(refDataInfo.prepareInfo, commonData, precursor, indicators, x, y, nx, ny, totalCells, optimizerValues, optimizerCondition); 
 %expectedOut=interface_expectedOutput(refDataInfo, commonData); 
 splitResult=0;
 nPoll=size(precursorNames, 1);
 intermediateResult=interface_prepare(refDataInfo, commonData); 
 % 20160418 : Quadrant version
 %if nPoll<>6 6;
 %if commonData.firstGuess
 %overwrite nPoll for FG 5 polls case
 %nPoll=6;
 expectedOut=interface_get_ExpectedOutput(intermediateResult, commonData, nx, ny, nPoll);
 % 20160418 : First Guess version
 %expectedOut=FG_interface_get_ExpectedOutput(intermediateResult, commonData, nx, ny, nPoll);
 %expectedOut=firstguess_get_expectedOutput(intermediateResult, commonData, nx, ny, nPoll);
 nextStartIndex=1;
 for i=1:nPoll
    precVals=precursorValues(:,i);
    %save('C:\data\work\projects\riat\precVal_wrong', 'precVals');
    % 20160418 : Mix versionQuadrant version
     [thisIntermediateResult, thisResult]=do_Job(refDataInfo, commonData, precVals, i, x, y, nx, ny, ...
     precursorNames(i), mathIndex, splitResult);
    % 20160418 : First Guess version
    %[thisIntermediateResult, thisResult]=firstguess_do_Job(refDataInfo, commonData, precVals, i, x, y, nx, ny, ...
    %precursorNames(i), mathIndex, splitResult);
    thisSize=size(thisResult.resGrid,1);
    expectedOut(nextStartIndex:nextStartIndex+thisSize-1)=thisResult.resGrid(1:thisSize);
%     imagesc(thisResult.resGrid(1:thisSize));
    nextStartIndex=nextStartIndex+thisSize;
    %par1_wrong=thisResult.resGrid;
    %save('C:\data\work\projects\riat\par1_wrong', 'par1_wrong');
 end
 %overall sparse
 %save('C:\data\work\projects\riat\all_bef_sparse_Wrong', 'expectedOut');
 if (doSparse==1) expectedOut=sparse(expectedOut); end
 %overall filter
 %save('C:\data\work\projects\riat\all_bef_filter_Wrong', 'expectedOut');
 if (doFilter==1) expectedOut=expectedOut((indoptrep==1 | indoptrep==2)); end
 %save('C:\data\work\projects\riat\all_final_Wrong', 'expectedOut');
 resultInfo.finalGrid=expectedOut;
 %resultInfo=interface_finalize(refDataInfo, intermediateResult, commonData, precursor, indicators, x, y, nx, ny, totalCells, optimizerValues, optimizerCondition);

end