function [ intermediateInfo  ] = interface_prepareData( refInfo, commonDataInfo, precursor, indicators, x, y, nx, ny, totalCells, ...
    optimizerValues, optimizerCondition)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% final version: from external file...

if (refInfo.isCoded)
stringToRun='';
stringToRun=strcat(stringToRun, refInfo.functionName);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'refInfo,');
stringToRun=strcat(stringToRun,'commonDataInfo,');
stringToRun=strcat(stringToRun,'precursor,');
stringToRun=strcat(stringToRun,'indicators,');
stringToRun=strcat(stringToRun,'x,');
stringToRun=strcat(stringToRun,'y,');
stringToRun=strcat(stringToRun,'nx,');
stringToRun=strcat(stringToRun,'ny,');
stringToRun=strcat(stringToRun,'totalCells,');
stringToRun=strcat(stringToRun,'optimizerValues,');
stringToRun=strcat(stringToRun,'optimizerCondition');
stringToRun=strcat(stringToRun,')');
[intermediateInfo]=eval(stringToRun);
else
for k=1:refInfo.linesNo
    %remember to assign resGrid somewhere in the lines
    % use evalC to avoid annoying output
    T=evalc(cell2mat(refInfo.codeLines(k)));
    %resGrid=resGrid(optimizerCondition);
end

end

