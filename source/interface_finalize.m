function [ finalizeInfo  ] = interface_finalize( refInfo, intermediateResult, commonData, precursor, indicators, x, y, nx, ny, ...
    precName, mathIndex, splitResult)
% function [ finalizeInfo  ] = interface_finalize( refInfo, intermediateResult, commonData, precursor, indicators, x, y, nx, ny, totalCells, ...
%     optimizerValues, optimizerCondition)

internalRefInfo=refInfo.finalizeInfo;

if (internalRefInfo.isCoded)
stringToRun='';
stringToRun=strcat(stringToRun, internalRefInfo.functionName);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'refInfo,');
stringToRun=strcat(stringToRun,'intermediateResult,');
stringToRun=strcat(stringToRun,'commonData,');
stringToRun=strcat(stringToRun,'precursor,');
stringToRun=strcat(stringToRun,'indicators,');
stringToRun=strcat(stringToRun,'x,');
stringToRun=strcat(stringToRun,'y,');
stringToRun=strcat(stringToRun,'nx,');
stringToRun=strcat(stringToRun,'ny,');
%stringToRun=strcat(stringToRun,'totalCells,');
%stringToRun=strcat(stringToRun,'optimizerValues,');
%stringToRun=strcat(stringToRun,'optimizerCondition');
stringToRun=strcat(stringToRun,'precName,');
stringToRun=strcat(stringToRun,'mathIndex, ');
stringToRun=strcat(stringToRun,'splitResult');
stringToRun=strcat(stringToRun,')');
[finalizeInfo]=eval(stringToRun);
else
for k=1:internalRefInfo.linesNo
    %remember to assign resGrid somewhere in the lines
    % use evalC to avoid annoying output
    T=evalc(cell2mat(internalRefInfo.codeLines(k)));
    %resGrid=resGrid(optimizerCondition);
end

end

